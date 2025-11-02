# OpenFOAM 运行时选择机制详解教程

## 目录
1. [概述](#概述)
2. [TypeName 和 Debug 宏的选择](#typename-和-debug-宏的选择)
3. [运行时选择表的工作原理](#运行时选择表的工作原理)
4. [完整示例：FastChemistryModel](#完整示例fastchemistrymodel)
5. [常见模式总结](#常见模式总结)

---

## 概述

OpenFOAM 的运行时选择机制（Runtime Selection Mechanism）允许程序在运行时根据配置文件动态选择和创建对象，而不需要在编译时确定具体的类型。这是一种工厂模式（Factory Pattern）的实现。

### 核心组件

1. **TypeName 宏** - 为类声明类型名和调试信息
2. **defineTypeNameAndDebug 宏** - 定义静态成员变量
3. **运行时选择表** - 存储构造函数指针的哈希表
4. **addToRunTimeSelectionTable 宏** - 向表中添加派生类

---

## TypeName 和 Debug 宏的选择

### 1. 基础：声明 vs 定义

在 C++ 中，静态成员变量需要**声明**和**定义**两个步骤：

```cpp
// 声明（在 .H 文件中）
class MyClass {
    static const word typeName;  // 仅声明
    static int debug;            // 仅声明
};

// 定义（在 .C 文件中）
const word MyClass::typeName("MyClass");  // 分配内存
int MyClass::debug(0);                    // 分配内存
```

### 2. TypeName 宏（仅声明）

位置：`OpenFOAM-10/src/OpenFOAM/db/typeInfo/className.H:65-67`

```cpp
// 宏定义
#define ClassName(TypeNameString)                \
    ClassNameNoDebug(TypeNameString);            \
    static int debug

// 展开 ClassNameNoDebug
#define ClassNameNoDebug(TypeNameString)         \
    static const char* typeName_() { return TypeNameString; }  \
    static const ::Foam::word typeName
```

**使用示例：**
```cpp
class basicFastChemistryModel {
public:
    TypeName("basicFastChemistryModel");
};
```

**展开后：**
```cpp
class basicFastChemistryModel {
public:
    // 内联函数，完整定义
    static const char* typeName_() { return "basicFastChemistryModel"; }

    // 仅声明！需要在 .C 文件中定义
    static const ::Foam::word typeName;
    static int debug;
};
```

⚠️ **注意**：`TypeName` 只提供声明，**不分配内存**！

### 3. defineTypeNameAndDebug 宏（定义，用于普通类）

位置：`OpenFOAM-10/src/OpenFOAM/db/typeInfo/className.H:119-121`

```cpp
// 宏定义
#define defineTypeNameAndDebug(Type, DebugSwitch)  \
    defineTypeName(Type);                          \
    defineDebugSwitch(Type, DebugSwitch)

// 进一步展开
#define defineTypeName(Type)                       \
    const ::Foam::word Type::typeName(Type::typeName_())
```

**使用示例：**
```cpp
// basicFastChemistryModel.C
namespace Foam
{
    defineTypeNameAndDebug(basicFastChemistryModel, 0);
}
```

**展开后：**
```cpp
namespace Foam
{
    // 定义 typeName 静态成员，分配内存
    const ::Foam::word basicFastChemistryModel::typeName(
        basicFastChemistryModel::typeName_()  // 调用内联函数返回 "basicFastChemistryModel"
    );

    // 定义 debug 静态成员，分配内存
    int basicFastChemistryModel::debug(0);
}
```

✅ **用途**：为**非模板类**定义 `typeName` 和 `debug` 静态成员变量。

### 4. defineTemplateTypeNameAndDebugWithName 宏（定义，用于模板类）

位置：`OpenFOAM-10/src/OpenFOAM/db/typeInfo/className.H:129-131`

```cpp
// 宏定义
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)  \
    defineTemplateTypeNameWithName(Type, Name);                          \
    defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)

// 进一步展开
#define defineTemplateTypeNameWithName(Type, Name)  \
    template<>                                      \
    const ::Foam::word Type::typeName(Name)
```

**使用示例（来自我们的代码）：**
```cpp
// odeChemistrySolvers.C
typedef OptSeulex<FastChemistryModel> OptSeulexFastChemistryModel;

defineTemplateTypeNameAndDebugWithName(
    OptSeulexFastChemistryModel,
    "OptSeulex<FastChemistryModel>",  // 自定义名称
    0
);
```

**展开后：**
```cpp
// 模板特化定义
template<>
const ::Foam::word OptSeulexFastChemistryModel::typeName(
    "OptSeulex<FastChemistryModel>"
);

template<>
int OptSeulexFastChemistryModel::debug(0);
```

✅ **用途**：为**模板类的特化版本**定义 `typeName` 和 `debug`，并可以**自定义类型名**。

### 5. 如何选择？

| 情况 | 使用宏 | 原因 |
|------|--------|------|
| 普通非模板类 | `defineTypeNameAndDebug(MyClass, 0)` | 直接定义静态成员 |
| 模板类的实例化 | `defineTemplateTypeNameAndDebugWithName(MyClass<T>, "MyClass<T>", 0)` | 需要模板特化语法 |
| 需要自定义类型名 | `defineTemplateTypeNameAndDebugWithName(...)` | 可以指定显示名称 |

### 6. 为什么在 makeChemistrySolver 中使用 defineTemplateTypeNameAndDebugWithName？

来自 `odeChemistrySolvers.C:67-83`：

```cpp
#define makeChemistrySolver(Solver)                                    \
    typedef Solver<FastChemistryModel> Solver##FastChemistryModel;     \
    defineTemplateTypeNameAndDebugWithName                             \
    (                                                                  \
        Solver##FastChemistryModel,                                    \
        (                                                              \
            word(Solver##FastChemistryModel::typeName_())              \
          + "<FastChemistryModel>"                                     \
        ).c_str(),                                                     \
        0                                                              \
    );                                                                 \
    addToRunTimeSelectionTable                                         \
    (                                                                  \
        basicFastChemistryModel,                                       \
        Solver##FastChemistryModel,                                    \
        thermo                                                         \
    )
```

**调用示例：**
```cpp
makeChemistrySolver(OptSeulex);
```

**展开后：**
```cpp
// 1. 创建 typedef
typedef OptSeulex<FastChemistryModel> OptSeulexFastChemistryModel;

// 2. 定义类型名和调试信息（模板特化）
template<>
const ::Foam::word OptSeulexFastChemistryModel::typeName(
    (word("OptSeulex") + "<FastChemistryModel>").c_str()
);
template<>
int OptSeulexFastChemistryModel::debug(0);

// 3. 添加到运行时选择表（后面详解）
basicFastChemistryModel::addthermoConstructorToTable<OptSeulexFastChemistryModel>
    addOptSeulexFastChemistryModelthermoConstructorTobasicFastChemistryModelTable_;
```

**为什么使用 defineTemplateTypeNameAndDebugWithName？**

1. **OptSeulex<FastChemistryModel> 是模板类实例化**
   - `OptSeulex` 本身是模板类
   - `OptSeulex<FastChemistryModel>` 是特化实例

2. **需要自定义类型名**
   - 我们希望类型名显示为 `"OptSeulex<FastChemistryModel>"`
   - 而不是 typedef 的名字 `"OptSeulexFastChemistryModel"`

3. **模板特化语法**
   - `template<>` 前缀是必需的，只有 `defineTemplateTypeNameAndDebugWithName` 提供此语法

---

## 运行时选择表的工作原理

### 1. 概念

运行时选择表本质上是一个 **HashTable**，存储了：
- **键（Key）**：类型名字符串（如 `"OptSeulex<FastChemistryModel>"`）
- **值（Value）**：构造函数指针（函数指针，用于创建对象）

### 2. 在基类中声明表

位置：`basicFastChemistryModel.H:109-116`

```cpp
class basicFastChemistryModel
{
public:
    //- Runtime type information
    TypeName("basicFastChemistryModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,                          // 返回类型
        basicFastChemistryModel,          // 基类名
        thermo,                           // 表的名称（构造函数签名标识）
        (const fluidReactionThermo& thermo),  // 构造函数参数列表（声明）
        (thermo)                          // 构造函数参数列表（调用）
    );
};
```

#### declareRunTimeSelectionTable 宏展开

位置：`OpenFOAM-10/src/OpenFOAM/db/runTimeSelection/construction/runTimeSelectionTables.H:46-129`

**展开后（简化版）：**
```cpp
class basicFastChemistryModel
{
public:
    // 1. 定义构造函数指针类型
    typedef autoPtr<basicFastChemistryModel> (*thermoConstructorPtr)
        (const fluidReactionThermo& thermo);

    // 2. 定义 HashTable 类型
    typedef HashTable<thermoConstructorPtr, word, string::hash>
        thermoConstructorTable;

    // 3. 声明指向 HashTable 的静态指针
    static thermoConstructorTable* thermoConstructorTablePtr_;

    // 4. 声明表的构造和析构函数
    static void constructthermoConstructorTables();
    static void destroythermoConstructorTables();

    // 5. 定义辅助类模板，用于向表中添加条目
    template<class basicFastChemistryModelType>
    class addthermoConstructorToTable
    {
    public:
        // 创建对象的静态函数
        static autoPtr<basicFastChemistryModel> New(const fluidReactionThermo& thermo)
        {
            return autoPtr<basicFastChemistryModel>(
                new basicFastChemistryModelType(thermo)
            );
        }

        // 构造函数：在表中注册
        addthermoConstructorToTable(
            const word& lookup = basicFastChemistryModelType::typeName
        )
        {
            // 确保表已创建
            constructthermoConstructorTables();

            // 插入到表中
            if (!thermoConstructorTablePtr_->insert(lookup, New))
            {
                std::cerr << "Duplicate entry " << lookup << std::endl;
            }
        }

        // 析构函数
        ~addthermoConstructorToTable()
        {
            destroythermoConstructorTables();
        }
    };
};
```

**关键点**：
- `addthermoConstructorToTable` 是一个**模板类**
- 它的**构造函数**会自动向表中添加条目
- 当创建这个类的全局对象时，构造函数在 `main()` 之前执行（全局对象初始化）

### 3. 在基类 .C 文件中定义表

位置：`basicFastChemistryModel.C:51`

```cpp
namespace Foam
{
    defineTypeNameAndDebug(basicFastChemistryModel, 0);
    defineRunTimeSelectionTable(basicFastChemistryModel, thermo);
}
```

#### defineRunTimeSelectionTable 宏展开

位置：`OpenFOAM-10/src/OpenFOAM/db/runTimeSelection/construction/runTimeSelectionTables.H:273-277`

```cpp
#define defineRunTimeSelectionTable(baseType, argNames)  \
    defineRunTimeSelectionTablePtr(baseType, argNames);  \
    defineRunTimeSelectionTableConstructor(baseType, argNames);  \
    defineRunTimeSelectionTableDestructor(baseType, argNames)
```

**展开后：**
```cpp
// 1. 定义静态指针（初始化为 nullptr）
basicFastChemistryModel::thermoConstructorTable*
    basicFastChemistryModel::thermoConstructorTablePtr_ = nullptr;

// 2. 定义表的构造函数
void basicFastChemistryModel::constructthermoConstructorTables()
{
    static bool constructed = false;
    if (!constructed)
    {
        constructed = true;
        basicFastChemistryModel::thermoConstructorTablePtr_
            = new basicFastChemistryModel::thermoConstructorTable;
    }
}

// 3. 定义表的析构函数
void basicFastChemistryModel::destroythermoConstructorTables()
{
    if (basicFastChemistryModel::thermoConstructorTablePtr_)
    {
        delete basicFastChemistryModel::thermoConstructorTablePtr_;
        basicFastChemistryModel::thermoConstructorTablePtr_ = nullptr;
    }
}
```

**作用**：
- 为静态 HashTable 指针分配内存
- 提供延迟初始化机制（第一次使用时创建）

### 4. 在派生类中添加到表

位置：`odeChemistrySolvers.C:78-83`

```cpp
addToRunTimeSelectionTable
(
    basicFastChemistryModel,           // 基类
    OptSeulexFastChemistryModel,       // 派生类（typedef）
    thermo                             // 构造函数签名标识
)
```

#### addToRunTimeSelectionTable 宏展开

位置：`OpenFOAM-10/src/OpenFOAM/db/runTimeSelection/construction/addToRunTimeSelectionTable.H:35-40`

```cpp
#define addToRunTimeSelectionTable(baseType,thisType,argNames)  \
    baseType::add##argNames##ConstructorToTable<thisType>       \
        add##thisType##argNames##ConstructorTo##baseType##Table_
```

**展开后：**
```cpp
// 创建全局对象（在 main() 之前初始化）
basicFastChemistryModel::addthermoConstructorToTable<OptSeulexFastChemistryModel>
    addOptSeulexFastChemistryModelthermoConstructorTobasicFastChemistryModelTable_;
```

**执行流程：**

1. **程序启动前**：创建全局对象
2. **全局对象构造**：调用 `addthermoConstructorToTable` 构造函数
3. **构造函数内部**：
   ```cpp
   addthermoConstructorToTable(
       const word& lookup = OptSeulexFastChemistryModel::typeName
   )
   {
       // 1. 确保表已创建
       constructthermoConstructorTables();

       // 2. 获取类型名
       // lookup = "OptSeulex<FastChemistryModel>"

       // 3. 插入到表中
       thermoConstructorTablePtr_->insert(
           "OptSeulex<FastChemistryModel>",  // 键
           &New                              // 值：指向 New() 函数的指针
       );
   }
   ```

4. **表中现在包含**：
   ```
   HashTable {
       "OptRodas34<FastChemistryModel>" → 函数指针（创建 OptRodas34<FastChemistryModel>）
       "OptRosenbrock34<FastChemistryModel>" → 函数指针（创建 OptRosenbrock34<FastChemistryModel>）
       "OptSeulex<FastChemistryModel>" → 函数指针（创建 OptSeulex<FastChemistryModel>）
   }
   ```

### 5. 运行时查找和创建对象

位置：`basicFastChemistryModelNew.C:84-115`（简化）

```cpp
autoPtr<basicFastChemistryModel> basicFastChemistryModel::New(
    const fluidReactionThermo& thermo
)
{
    // 1. 读取配置文件
    IOdictionary chemistryDict(...);

    // 2. 获取求解器名称（如 "OptSeulex"）
    const word solverName = chemistryDict.lookup("solver");

    // 3. 获取方法名称（如 "FastChemistryModel"）
    const word methodName = chemistryDict.lookup("method");

    // 4. 构造查找键
    const word chemSolverNameName =
        solverName + '<' + methodName + '>';
    // 结果："OptSeulex<FastChemistryModel>"

    // 5. 在表中查找
    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(chemSolverNameName);

    // 6. 调用构造函数指针
    return autoPtr<basicFastChemistryModel>(
        cstrIter()(thermo)  // 调用找到的函数指针
    );
}
```

**查找流程：**
```
chemistryProperties 文件:
{
    solver   OptSeulex;
    method   FastChemistryModel;
}
       ↓
构造键: "OptSeulex<FastChemistryModel>"
       ↓
在 HashTable 中查找
       ↓
找到函数指针
       ↓
调用函数指针创建对象
       ↓
返回 autoPtr<OptSeulex<FastChemistryModel>>
（但类型擦除为 autoPtr<basicFastChemistryModel>）
```

---

## 完整示例：FastChemistryModel

### 文件结构

```
src/
├── basicFastChemistryModel/
│   ├── basicFastChemistryModel.H       # 基类声明，声明运行时选择表
│   ├── basicFastChemistryModel.C       # 基类实现，定义运行时选择表
│   └── basicFastChemistryModelNew.C    # 工厂方法，查表创建对象
├── ChemistryModel/
│   ├── FastChemistryModel.H            # 派生类声明
│   └── FastChemistryModel.C            # 派生类实现
├── ODE/
│   ├── OptSeulex.H                     # ODE 求解器（模板类）
│   └── OptSeulex.C
└── odeChemistrySolvers.C               # 注册所有求解器到表
```

### 第 1 步：基类声明运行时选择表

**文件：`basicFastChemistryModel.H`**

```cpp
namespace Foam
{

class basicFastChemistryModel
:
    public IOdictionary
{
public:
    //- Runtime type information
    TypeName("basicFastChemistryModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicFastChemistryModel,
        thermo,
        (const fluidReactionThermo& thermo),
        (thermo)
    );

    // Selectors
    static autoPtr<basicFastChemistryModel> New(
        const fluidReactionThermo& thermo
    );

    //- Destructor
    virtual ~basicFastChemistryModel();
};

}
```

**关键点：**
- `TypeName("basicFastChemistryModel")` 声明类型名
- `declareRunTimeSelectionTable(...)` 声明运行时选择表
- `New()` 是工厂方法

### 第 2 步：基类定义运行时选择表

**文件：`basicFastChemistryModel.C`**

```cpp
namespace Foam
{
    // 定义类型名和调试信息
    defineTypeNameAndDebug(basicFastChemistryModel, 0);

    // 定义运行时选择表（创建 HashTable）
    defineRunTimeSelectionTable(basicFastChemistryModel, thermo);
}

basicFastChemistryModel::basicFastChemistryModel(
    const fluidReactionThermo& thermo
)
:
    IOdictionary(...),
    mesh_(thermo.T().mesh()),
    thermo_(thermo)
{}

basicFastChemistryModel::~basicFastChemistryModel()
{}
```

**关键点：**
- `defineTypeNameAndDebug` 定义 `typeName` 和 `debug` 静态成员
- `defineRunTimeSelectionTable` 创建 HashTable 及相关函数

### 第 3 步：派生类声明

**文件：`FastChemistryModel.H`**

```cpp
namespace Foam
{

class FastChemistryModel
:
    public basicFastChemistryModel
{
public:
    //- Runtime type information
    TypeName("FastChemistryModel");

    // Constructors
    FastChemistryModel(const fluidReactionThermo& thermo);

    //- Destructor
    virtual ~FastChemistryModel();
};

}
```

**关键点：**
- 继承自 `basicFastChemistryModel`
- 声明自己的类型名

### 第 4 步：派生类定义

**文件：`FastChemistryModel.C`**

```cpp
namespace Foam
{
    // 定义类型名和调试信息
    defineTypeNameAndDebug(FastChemistryModel, 0);
}

FastChemistryModel::FastChemistryModel(
    const fluidReactionThermo& thermo
)
:
    basicFastChemistryModel(thermo),
    reaction(...)
{
    // 初始化
}

FastChemistryModel::~FastChemistryModel()
{}
```

### 第 5 步：注册模板类到运行时选择表

**文件：`odeChemistrySolvers.C`**

```cpp
#include "OptRosenbrock34.H"
#include "OptRodas34.H"
#include "OptSeulex.H"
#include "FastChemistryModel.H"
#include "basicFastChemistryModel.H"

namespace Foam
{

// 定义 FastChemistryModel 的类型名和调试信息
defineTypeNameAndDebug(FastChemistryModel, 0);

// 定义宏，用于注册 ODE 求解器
#define makeChemistrySolver(Solver)                                    \
    /* 创建 typedef */                                                 \
    typedef Solver<FastChemistryModel> Solver##FastChemistryModel;     \
    /* 定义模板类的类型名和调试信息 */                                  \
    defineTemplateTypeNameAndDebugWithName                             \
    (                                                                  \
        Solver##FastChemistryModel,                                    \
        (                                                              \
            word(Solver##FastChemistryModel::typeName_())              \
          + "<FastChemistryModel>"                                     \
        ).c_str(),                                                     \
        0                                                              \
    );                                                                 \
    /* 添加到运行时选择表 */                                           \
    addToRunTimeSelectionTable                                         \
    (                                                                  \
        basicFastChemistryModel,                                       \
        Solver##FastChemistryModel,                                    \
        thermo                                                         \
    )

// 注册所有 ODE 求解器
makeChemistrySolver(OptRosenbrock34);
makeChemistrySolver(OptRodas34);
makeChemistrySolver(OptSeulex);

}
```

**展开 `makeChemistrySolver(OptSeulex)` 后：**

```cpp
// 1. 创建 typedef
typedef OptSeulex<FastChemistryModel> OptSeulexFastChemistryModel;

// 2. 定义类型名（模板特化）
template<>
const ::Foam::word OptSeulexFastChemistryModel::typeName(
    "OptSeulex<FastChemistryModel>"
);

template<>
int OptSeulexFastChemistryModel::debug(0);

// 3. 添加到运行时选择表（创建全局对象）
basicFastChemistryModel::addthermoConstructorToTable<OptSeulexFastChemistryModel>
    addOptSeulexFastChemistryModelthermoConstructorTobasicFastChemistryModelTable_;
```

### 第 6 步：工厂方法（查表创建对象）

**文件：`basicFastChemistryModelNew.C`**

```cpp
Foam::autoPtr<Foam::basicFastChemistryModel>
Foam::basicFastChemistryModel::New(
    const fluidReactionThermo& thermo
)
{
    // 读取配置字典
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // 获取求解器名称
    const word solverName = chemistryDict.lookup("solver");

    // 获取方法名称
    const word methodName = chemistryDict.lookup("method");

    // 构造查找键
    const word chemSolverNameName =
        solverName + '<' + methodName + '>';

    Info<< "Looking for: " << chemSolverNameName << nl;

    // 在表中查找
    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(chemSolverNameName);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        // 未找到，打印可用选项
        FatalErrorInFunction
            << "Unknown solver " << chemSolverNameName << nl << nl
            << "Valid solvers are:" << nl
            << thermoConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    // 调用构造函数指针创建对象
    return autoPtr<basicFastChemistryModel>(cstrIter()(thermo));
}
```

### 第 7 步：运行时配置

**文件：`constant/chemistryProperties`**

```cpp
chemistryType
{
    solver          OptSeulex;
    method          FastChemistryModel;
}

chemistry       on;
initialChemicalTimeStep 1e-8;
```

### 运行时执行流程

```
程序启动
    ↓
全局对象初始化（在 main() 之前）
    ↓
addOptSeulexFastChemistryModelthermoConstructorTobasicFastChemistryModelTable_
构造函数执行
    ↓
向 thermoConstructorTablePtr_ 插入条目：
    键："OptSeulex<FastChemistryModel>"
    值：函数指针（创建 OptSeulex<FastChemistryModel> 对象）
    ↓
main() 执行
    ↓
调用 basicFastChemistryModel::New(thermo)
    ↓
读取 chemistryProperties，获得：
    solver = "OptSeulex"
    method = "FastChemistryModel"
    ↓
构造查找键："OptSeulex<FastChemistryModel>"
    ↓
在 thermoConstructorTablePtr_ 中查找
    ↓
找到对应的函数指针
    ↓
调用函数指针：cstrIter()(thermo)
    ↓
创建 OptSeulex<FastChemistryModel> 对象
    ↓
返回 autoPtr<basicFastChemistryModel>
（实际指向 OptSeulex<FastChemistryModel> 对象）
```

---

## 常见模式总结

### 模式 1：普通非模板类

```cpp
// .H 文件
class MyClass : public BaseClass
{
public:
    TypeName("MyClass");
};

// .C 文件
defineTypeNameAndDebug(MyClass, 0);
addToRunTimeSelectionTable(BaseClass, MyClass, argNames);
```

### 模式 2：模板类的特化

```cpp
// .H 文件
template<class T>
class MyTemplateClass : public BaseClass
{
public:
    ClassName("MyTemplateClass");
};

// .C 文件或注册文件
typedef MyTemplateClass<ConcreteType> MyConcreteClass;

defineTemplateTypeNameAndDebugWithName(
    MyConcreteClass,
    "MyTemplateClass<ConcreteType>",
    0
);

addToRunTimeSelectionTable(BaseClass, MyConcreteClass, argNames);
```

### 模式 3：使用宏批量注册

```cpp
// 定义注册宏
#define makeMyClass(SolverType)                                \
    typedef SolverType<MyModel> SolverType##MyModel;           \
    defineTemplateTypeNameAndDebugWithName(                    \
        SolverType##MyModel,                                   \
        word(SolverType##MyModel::typeName_()) + "<MyModel>",  \
        0                                                      \
    );                                                         \
    addToRunTimeSelectionTable(BaseClass, SolverType##MyModel, argNames)

// 批量注册
makeMyClass(Solver1);
makeMyClass(Solver2);
makeMyClass(Solver3);
```

---

## 关键要点总结

### 1. 声明 vs 定义

| 组件 | 声明位置 | 定义位置 |
|------|----------|----------|
| `TypeName` 宏 | 类的 .H 文件 | - |
| `typeName` 静态成员 | `TypeName` 宏中 | `defineTypeNameAndDebug` 宏中 |
| `debug` 静态成员 | `TypeName` 宏中 | `defineTypeNameAndDebug` 宏中 |
| 运行时选择表 | `declareRunTimeSelectionTable` | `defineRunTimeSelectionTable` |

### 2. 普通类 vs 模板类

| 类型 | 使用宏 |
|------|--------|
| 普通非模板类 | `defineTypeNameAndDebug` |
| 模板类特化 | `defineTemplateTypeNameAndDebugWithName` |

### 3. 运行时选择表的三要素

1. **声明**：`declareRunTimeSelectionTable` 在基类 .H 中
2. **定义**：`defineRunTimeSelectionTable` 在基类 .C 中
3. **注册**：`addToRunTimeSelectionTable` 在派生类的注册文件中

### 4. 全局对象初始化顺序

```
程序加载
    ↓
全局对象构造（main() 之前）
    ↓
addToRunTimeSelectionTable 创建全局对象
    ↓
全局对象的构造函数向表中添加条目
    ↓
main() 开始执行
    ↓
表已经包含所有注册的类型
```

### 5. 为什么可以在运行时选择？

- **编译时**：所有派生类都注册到同一个全局 HashTable
- **运行时**：根据配置文件的字符串查找 HashTable
- **多态**：通过基类指针调用派生类的虚函数

---

## 调试技巧

### 查看注册的所有类型

在 `New()` 函数中添加：

```cpp
Info<< "Available types: " << nl
    << thermoConstructorTablePtr_->sortedToc() << endl;
```

### 检查类型名是否正确

使用 `nm` 命令查看符号：

```bash
nm libMyLibrary.so | grep typeName
```

应该看到：
```
B MyClass::typeName        # B = BSS (已定义)
U MyClass::typeName        # U = Undefined（未定义，链接错误）
```

### 检查表是否包含条目

在库加载后添加调试输出：

```cpp
if (thermoConstructorTablePtr_)
{
    Info<< "Table size: " << thermoConstructorTablePtr_->size() << endl;
    Info<< "Entries: " << thermoConstructorTablePtr_->toc() << endl;
}
```

---

## 参考资料

- OpenFOAM 源码：`src/OpenFOAM/db/typeInfo/className.H`
- OpenFOAM 源码：`src/OpenFOAM/db/runTimeSelection/construction/runTimeSelectionTables.H`
- OpenFOAM 源码：`src/OpenFOAM/db/runTimeSelection/construction/addToRunTimeSelectionTable.H`
- 本项目示例：`src/basicFastChemistryModel/`
- 本项目示例：`src/odeChemistrySolvers.C`

---

**教程结束**

希望这个教程帮助你理解了 OpenFOAM 运行时选择机制的工作原理！
