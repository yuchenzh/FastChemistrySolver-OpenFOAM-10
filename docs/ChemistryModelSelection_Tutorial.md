# OpenFOAM 化学模型选择机制教程

## 概述

本教程详细说明 OpenFOAM 中 `basicChemistryModel` 的运行时选择机制（Run-Time Selection Table, RTST）如何工作，以及如何正确配置和使用自定义化学模型。

## 1. 选择机制的核心原理

### 1.1 运行时选择表（RTST）

OpenFOAM 使用运行时选择表机制来动态选择化学模型，这允许用户在运行时通过配置文件选择不同的求解器和方法，而无需重新编译代码。

**关键文件：**
- `src/basicChemistryModel/basicChemistryModel.H:108-116` - 声明运行时选择表
- `src/basicChemistryModel/basicChemistryModel.C:51` - 定义运行时选择表
- `src/basicChemistryModel/basicChemistryModelNew.C` - 实现选择器函数

### 1.2 选择表的声明

在 `basicChemistryModel.H` 中：

```cpp
declareRunTimeSelectionTable
(
    autoPtr,
    basicChemistryModel,
    thermo,
    (const fluidReactionThermo& thermo),
    (thermo)
);
```

**参数说明：**
- `autoPtr` - 返回类型（智能指针）
- `basicChemistryModel` - 基类名称
- `thermo` - 选择表的名称（通常基于构造函数参数类型）
- `(const fluidReactionThermo& thermo)` - 构造函数参数
- `(thermo)` - 参数名称

### 1.3 选择表的定义

在 `basicChemistryModel.C` 中：

```cpp
namespace Foam
{
    defineTypeNameAndDebug(basicChemistryModel, 0);
    defineRunTimeSelectionTable(basicChemistryModel, thermo);
}
```

## 2. 选择器函数的工作流程

### 2.1 完整的选择过程

`basicChemistryModelNew.C` 中的 `New()` 函数实现了以下步骤：

#### 步骤 1: 读取配置文件

```cpp
IOdictionary chemistryDict
(
    IOobject
    (
        thermo.phasePropertyName("chemistryProperties"),  // 文件名
        thermo.T().mesh().time().constant(),               // 路径: constant/
        thermo.T().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);
```

**读取位置：** `constant/chemistryProperties`

#### 步骤 2: 检查配置格式

```cpp
if (!chemistryDict.isDict("chemistryType"))
{
    FatalErrorInFunction
        << "Template parameter based chemistry solver selection is no "
        << "longer supported. Please create a chemistryType dictionary"
        << "instead."
        << exit(FatalError);
}
```

**要求：** 必须使用 `chemistryType` 字典格式（不再支持旧的模板参数格式）

#### 步骤 3: 提取求解器和方法名称

```cpp
const dictionary& chemistryTypeDict = chemistryDict.subDict("chemistryType");

const word solverName = chemistryTypeDict.lookupBackwardsCompatible<word>
(
    {"solver", "chemistrySolver"}  // 支持两种关键字
);

const word methodName = chemistryTypeDict.lookupOrDefault<word>
(
    "method",
    "chemistryModel"  // 默认值
);
```

#### 步骤 4: 构造完整的类型名称

```cpp
const word chemSolverNameName =
    solverName + '<' + methodName + '<' + thermo.thermoName() + ">>";
```

**示例：**
- `solverName` = "OptRodas34"
- `methodName` = "FastChemistryModel"
- `thermo.thermoName()` = "sutherland<janaf<perfectGas>,sensibleInternalEnergy>"
- **结果：** `"OptRodas34<FastChemistryModel<sutherland<janaf<perfectGas>,sensibleInternalEnergy>>>"`

#### 步骤 5: 在选择表中查找

```cpp
typename thermoConstructorTable::iterator cstrIter =
    thermoConstructorTablePtr_->find(chemSolverNameName);

if (cstrIter == thermoConstructorTablePtr_->end())
{
    // 未找到：尝试动态编译或报错
}
```

#### 步骤 6: 动态编译（如果启用）

```cpp
if (dynamicCode::allowSystemOperations
 && !dynamicCode::resolveTemplate(basicChemistryModel::typeName).empty())
{
    List<Pair<word>> substitutions
    (
        basicThermo::thermoNameComponents(thermo.thermoName())
    );

    substitutions.append({"solver", solverName});
    substitutions.append({"method", methodName});

    compileTemplate chemistryModel
    (
        basicChemistryModel::typeName,
        chemSolverNameName,
        substitutions
    );

    cstrIter = thermoConstructorTablePtr_->find(chemSolverNameName);
}
```

**功能：** 如果预编译的类型不存在，系统会尝试从模板文件动态编译

#### 步骤 7: 创建对象实例

```cpp
return autoPtr<basicChemistryModel>(cstrIter()(thermo));
```

## 3. 配置文件详解

### 3.1 基本配置格式

**文件位置：** `constant/chemistryProperties`

```cpp
chemistryType
{
    solver          OptRodas34;          // ODE求解器名称
    method          FastChemistryModel;  // 化学模型方法
}

chemistry       on;                      // 启用/禁用化学反应
jacobian        exact;                   // Jacobian 类型：fast/exact

initialChemicalTimeStep  1e-8;          // 初始化学时间步长
maxChemicalTimeStep      great;         // 最大化学时间步长（可选）

// 求解器特定的系数
OptRodas34Coeffs
{
    absTol          1e-12;              // 绝对容差
    relTol          1e-6;               // 相对容差
}
```

### 3.2 常用配置选项

#### solver 选项：
- `ode` - 标准 ODE 求解器
- `OptRodas34` - 优化的 Rodas34 求解器
- `OptRosenbrock34` - 优化的 Rosenbrock34 求解器
- `OptSeulex` - 优化的 Seulex 求解器

#### method 选项：
- `chemistryModel` - 标准化学模型（默认）
- `FastChemistryModel` - 快速化学模型（本项目）

#### jacobian 选项：
- `fast` - 快速近似 Jacobian
- `exact` - 精确 Jacobian（推荐用于准确性）

## 4. 热力学模型名称解析

### 4.1 热力学名称的组成

热力学模型名称（`thermo.thermoName()`）由多个组件组成：

```
<transport><thermo><equationOfState><specie><energy>
```

**示例：**
```
sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>
```

**组件说明：**
1. `transport` - 传输模型（如 sutherland, const）
2. `thermo` - 热力学数据（如 janaf, hConst）
3. `equationOfState` - 状态方程（如 perfectGas, incompressiblePerfectGas）
4. `specie` - 物种模型
5. `energy` - 能量类型（如 sensibleInternalEnergy, sensibleEnthalpy）

### 4.2 完整类型名称示例

最终构造的类型名称：
```
OptRodas34<FastChemistryModel<sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>>>
```

这个名称会被用于在运行时选择表中查找对应的构造函数。

## 5. 如何添加自定义化学模型

### 5.1 创建派生类

```cpp
// FastChemistryModel.H
template<class ThermoType>
class FastChemistryModel
:
    public ode<chemistryModel<ThermoType>>
{
    // ... 类定义

public:

    //- Runtime type information
    TypeName("FastChemistryModel");

    // 构造函数
    FastChemistryModel(const fluidReactionThermo& thermo);
};
```

### 5.2 注册到选择表

在实现文件中（通常是 `.C` 文件）：

```cpp
// 为每种支持的热力学类型注册
#include "makeChemistrySolverType.H"

namespace Foam
{
    // 定义具体实例化
    typedef chemistryModel<rhoThermo::thermoType> rhoChemistryModel;
    typedef FastChemistryModel<rhoThermo::thermoType> FastRhoChemistryModel;

    // 添加到选择表
    addChemistrySolverToRunTimeSelectionTable
    (
        FastRhoChemistryModel,
        rhoChemistryModel
    );
}
```

### 5.3 编译系统配置

在 `Make/files` 中：

```bash
odeChemistrySolvers.C

LIB = $(FOAM_USER_LIBBIN)/libFastChemistryModel
```

在 `Make/options` 中：

```bash
EXE_INC = \
    -I$(LIB_SRC)/thermophysicalModels/reactionThermo/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/chemistryModel/lnInclude \
    -I$(LIB_SRC)/ODE/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

LIB_LIBS = \
    -lreactionThermophysicalModels \
    -lspecie \
    -lchemistryModel \
    -lODE \
    -lfiniteVolume \
    -lmeshTools
```

## 6. 调试和故障排查

### 6.1 常见错误

#### 错误 1: 类型未找到

```
Unknown basicChemistryModel type
{
    solver OptRodas34;
    method FastChemistryModel;
}
```

**原因：**
- 库未正确编译
- 库未在求解器中链接
- 类型名称拼写错误

**解决方法：**
1. 检查库是否编译：`ls $FOAM_USER_LIBBIN/libFastChemistryModel.so`
2. 在 `system/controlDict` 中添加库：
   ```cpp
   libs
   (
       "libFastChemistryModel.so"
   );
   ```

#### 错误 2: 声明与定义不匹配

```
error: no declaration matches 'Foam::scalarField Foam::FastChemistryModel<ThermoType>::getRRGivenYTP(...)'
```

**原因：**
- 头文件中函数声明为 `const`，但实现中缺少 `const`
- 参数类型不匹配

**解决方法：**
确保函数签名完全一致，包括 `const` 限定符

### 6.2 查看可用的类型

当发生类型未找到错误时，OpenFOAM 会列出所有可用的组合：

```
Valid solver/method combinations for this thermodynamic model are:

solver           method
ode              chemistryModel
OptRodas34       chemistryModel
OptRodas34       FastChemistryModel
```

## 7. 完整示例

### 7.1 配置文件示例

**constant/chemistryProperties:**

```cpp
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      chemistryProperties;
}

chemistryType
{
    solver          OptRodas34;
    method          FastChemistryModel;
}

chemistry       on;
jacobian        exact;
initialChemicalTimeStep 1e-8;

OptRodas34Coeffs
{
    absTol          1e-12;
    relTol          1e-6;
}

#include "reactions"
```

### 7.2 控制字典示例

**system/controlDict:**

```cpp
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     reactingFoam;

libs
(
    "libFastChemistryModel.so"  // 加载自定义库
);

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         0.001;
deltaT          1e-6;
writeControl    adjustable;
writeInterval   1e-5;
```

## 8. 高级主题

### 8.1 动态编译

OpenFOAM 支持动态编译未预编译的化学模型类型。这需要：

1. 模板文件存在于正确位置
2. 系统允许运行时编译（`allowSystemOperations = true`）
3. 提供正确的替换参数

### 8.2 模板实例化

对于模板类，需要为每种支持的热力学类型显式实例化：

```cpp
// 在 odeChemistrySolvers.C 中
#include "makeChemistrySolverType.H"
#include "thermoPhysicsTypes.H"
#include "FastChemistryModel.H"

// 为所有标准热力学类型创建实例
makeChemistrySolverTypes(FastChemistryModel, psiChemistryModel);
makeChemistrySolverTypes(FastChemistryModel, rhoChemistryModel);
```

### 8.3 多相流系统

对于多相流，配置文件名称会包含相名：

```cpp
// 对于相 "gas"
constant/gas.chemistryProperties

// 或者使用 phaseProperties
constant/phaseProperties
```

## 9. 性能考虑

### 9.1 选择合适的求解器

- **ode**: 通用，但可能较慢
- **OptRodas34**: 对刚性系统高效
- **OptSeulex**: 对非刚性系统高效

### 9.2 Jacobian 选择

- **exact**: 更准确，但计算成本高
- **fast**: 近似，但计算快速

### 9.3 容差设置

```cpp
absTol  1e-12;  // 绝对容差：更小 = 更准确但更慢
relTol  1e-6;   // 相对容差：更小 = 更准确但更慢
```

## 10. 总结

OpenFOAM 的化学模型选择机制提供了：

1. **灵活性**：运行时选择不同的求解器和方法
2. **可扩展性**：轻松添加新的化学模型
3. **类型安全**：编译时类型检查
4. **动态编译**：支持运行时模板实例化

**关键要点：**
- 配置文件位于 `constant/chemistryProperties`
- 使用 `chemistryType` 字典指定求解器和方法
- 完整类型名称由求解器、方法和热力学模型组合而成
- 自定义模型需要正确注册到运行时选择表
- 确保库在 `system/controlDict` 中被加载

## 参考资源

- OpenFOAM 用户指南：https://openfoam.org/resources/
- OpenFOAM 程序员指南
- 源代码：`$FOAM_SRC/thermophysicalModels/chemistryModel/`
