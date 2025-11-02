# ChemistryModelHelper 快速参考

## 概述

`ChemistryModelHelper` 是一个轻量级工具类，允许你调用 `FastChemistryModel` 的扩展功能（如 `getRRGivenYTP`），而无需修改 `basicChemistryModel` 基类。

**特点**：
- ✅ 极简设计：只转换一次，后续调用零开销
- ✅ 高性能：使用 `static_cast` 直接转换
- ✅ 不修改基类：保持代码模块化

## 快速开始

### 1. 在求解器中使用

```cpp
#include "ChemistryModelHelper.H"

// 创建化学模型
autoPtr<basicChemistryModel> chemistry = basicChemistryModel::New(thermo);

// 创建 helper（类型转换只执行一次）
ChemistryModelHelper helper(chemistry());

// 使用 helper 调用 getRRGivenYTP（零开销）
scalarField RR = helper.getRRGivenYTP
(
    Y,              // 物种质量分数
    T,              // 温度 [K]
    p,              // 压力 [Pa]
    dt,             // 时间步长 [s]
    deltaTChem,     // 化学时间步长（输出）[s]
    rho,            // 密度 [kg/m^3]
    rho             // 参考密度 [kg/m^3]
);
```

### 2. 编译配置

**Make/options:**

```makefile
EXE_INC = \
    -I../../src/ChemistryModel \
    -I../../src/lnInclude \
    ...

EXE_LIBS = \
    -L$(FOAM_USER_LIBBIN) \
    -lFastChemistryModel \
    ...
```

### 3. 配置文件

**constant/chemistryProperties:**

```cpp
chemistryType
{
    solver          OptRodas34;
    method          FastChemistryModel;  // 必须使用 FastChemistryModel
}
```

## 编译和运行

```bash
# 1. 加载 OpenFOAM 环境
source /home/zhouyuchen/OpenFOAM/of/OpenFOAM-10/OpenFOAM-10/etc/bashrc

# 2. 编译库
cd src
wclean && wmake -j

# 3. 编译求解器
cd ../testSolver
wclean && wmake -j

# 4. 运行
cd ../Sandia
testFoam
```

## 文件位置

- **项目记忆**: `.claude/CLAUDE.md`
- **Helper 类**: `src/ChemistryModel/ChemistryModelHelper.H` 和 `ChemistryModelHelperI.H`
- **详细教程**: `docs/ChemistryModelHelper_Tutorial.md`
- **选择机制**: `docs/ChemistryModelSelection_Tutorial.md`

## 常见问题

### Q: 编译时找不到 ChemistryModelHelper.H？

A: 检查 `Make/options` 中的包含路径：
```makefile
-I../../src/ChemistryModel
-I../../src/lnInclude
```

### Q: 链接错误 undefined reference？

A: 确保 `Make/options` 中链接了库：
```makefile
-L$(FOAM_USER_LIBBIN) -lFastChemistryModel
```

### Q: 运行时崩溃或结果错误？

A: 确保 `constant/chemistryProperties` 中使用 FastChemistryModel：
```cpp
chemistryType
{
    method  FastChemistryModel;  // 必须！
}
```

## 性能特点

### 设计原则
- **一次转换**：在 Helper 构造时执行类型转换并缓存指针
- **直接调用**：使用 `static_cast`，无运行时类型检查
- **零开销抽象**：后续调用 `getRRGivenYTP` 接近零开销

### 性能对比

| 调用方式 | 开销 | 说明 |
|---------|------|------|
| 直接调用（基类无方法） | N/A | ❌ 不可用 |
| dynamic_cast 每次检查 | 中等 | RTTI 开销 |
| **Helper 类（本方案）** | **≈0** | ✅ 推荐 |

### 最佳实践

如果需要多次调用，在循环外创建 helper：

```cpp
// 在循环外创建（只转换一次）
ChemistryModelHelper helper(chemistry());

while (runTime.loop())
{
    // 多次调用，无重复转换开销
    scalarField RR = helper.getRRGivenYTP(...);
}
```

## 技术细节

### 支持的热力学类型
FastChemistryModel 只支持一种类型：
```cpp
sutherlandTransport<sensibleEnthalpy<janafThermo<perfectGas<specie>>>>
```

见 `src/odeChemistrySolvers.C:63`

### 实现方式
1. **构造函数**：执行类型转换并缓存指针
2. **getRRGivenYTP**：直接使用缓存的指针调用

详见 `docs/ChemistryModelHelper_Simple.md`

## 更多信息

- `docs/ChemistryModelHelper_Simple.md` - 简化版实现说明
- `docs/ChemistryModelSelection_Tutorial.md` - 化学模型选择机制
- `.claude/CLAUDE.md` - 项目配置信息
