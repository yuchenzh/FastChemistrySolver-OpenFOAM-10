# FastChemistrySolver-OpenFOAM-10 项目信息

## 项目概述
这是一个基于 OpenFOAM-10 的快速化学求解器库，提供优化的化学反应求解功能。

## 环境配置

### OpenFOAM 环境加载
```bash
source /home/zhouyuchen/OpenFOAM/of/OpenFOAM-10/OpenFOAM-10/etc/bashrc
```

### 编译命令
```bash
# 清理并编译库
wclean && wmake -j

# 清理并编译求解器
cd testSolver
wclean && wmake -j
```

### 运行测试
```bash
cd Sandia  # 或其他算例目录
testFoam   # 运行测试求解器
```

## 项目结构

### 核心源代码
- `src/ChemistryModel/` - FastChemistryModel 核心实现
  - `FastChemistryModel.H` - 头文件
  - `FastChemistryModel.C` - 实现文件

- `src/basicChemistryModel/` - 化学模型基类
  - `basicChemistryModel.H` - 基类头文件
  - `basicChemistryModel.C` - 基类实现
  - `basicChemistryModelNew.C` - 运行时选择机制

- `src/ODE/` - ODE 求解器
  - `OptRodas34/` - 优化的 Rodas34 求解器
  - `OptRosenbrock34/` - 优化的 Rosenbrock34 求解器
  - `OptSeulex/` - 优化的 Seulex 求解器

### 测试求解器
- `testSolver/` - 测试求解器目录
  - `reactingFoam.C` - 求解器主文件
  - `Make/files` - 编译文件列表
  - `Make/options` - 编译选项

### 算例目录
- `tutorial/Ignition0D/` - 0维点火算例
- `tutorial/CounterFlowFlame2D/` - 2维对撞火焰算例
- `tutorial/LES/` - LES算例
- `tutorial/RAS/` - RAS算例
- `Sandia/` - Sandia火焰算例

## 关键特性

### FastChemistryModel 特性
- 快速化学反应率计算
- 支持多种 ODE 求解器
- 负载均衡功能
- 优化的 Jacobian 计算（fast/exact）

### 配置文件
位置：`constant/chemistryProperties`

```cpp
chemistryType
{
    solver          OptRodas34;         // ODE求解器
    method          FastChemistryModel; // 化学模型方法
}

chemistry       on;
jacobian        exact;  // fast 或 exact
initialChemicalTimeStep 1e-8;
```

## 当前任务

### 总体目标
将 FastChemistrySolver 集成到 OpenFOAM 标准化学模型框架，使其能够在其他求解器和库中使用。

### 实施步骤（按顺序执行）

#### 步骤 1: 重命名 basicChemistryModel ⏳
**目标**: 避免与 OpenFOAM 标准库的 `basicChemistryModel` 冲突

**任务**:
- 将自定义的 `basicChemistryModel` 重命名为新名称（如 `fastBasicChemistryModel`）
- 更新所有相关文件中的引用
- 修复所有依赖关系
- 确保编译通过

**涉及文件**:
- `src/basicChemistryModel/basicChemistryModel.{H,C}`
- `src/basicChemistryModel/basicChemistryModelNew.C`
- `src/ChemistryModel/FastChemistryModel.{H,C}`
- 所有包含该类的头文件

#### 步骤 2: 修改 chemistryModelNew 构造逻辑 ⏳
**目标**: 实现精确构造，从子字典读取配置而非整个 chemistryProperties

**当前问题**:
- 构造函数读取整个 `chemistryProperties` 字典
- 需要改为读取其中的特定子字典

**修改方案**:
- 修改 `basicChemistryModelNew::New()` 函数
- 从 `chemistryProperties` 的某个子字典读取配置
- 支持更精确的模型选择和参数配置

**涉及文件**:
- `src/basicChemistryModel/basicChemistryModelNew.C`

#### 步骤 3: 添加 getRRGivenYTP 到基类 ⏳
**目标**: 将 `getRRGivenYTP` 作为基类虚函数，移除临时的 Helper 类

**任务**:
- 在重命名后的基类中添加虚函数 `getRRGivenYTP`
- `FastChemistryModel` 重写该方法
- 移除 `ChemistryModelHelper` 类
- 更新所有调用代码直接使用基类接口

**涉及文件**:
- 基类头文件（重命名后）
- `src/ChemistryModel/FastChemistryModel.{H,C}`
- `src/ChemistryModel/ChemistryModelHelper.{H}` (删除)
- `src/ChemistryModel/ChemistryModelHelperI.H` (删除)

### 当前进度
- ✅ ChemistryModelHelper 临时实现完成（临时方案）
- ⏳ **进行中**: 步骤 1 - 重命名 basicChemistryModel
- ⏳ **待完成**: 步骤 2 - 修改构造逻辑
- ⏳ **待完成**: 步骤 3 - 添加 getRRGivenYTP 到基类

### 设计原则
1. 遵循 OpenFOAM 运行时选择表 (RTST) 机制
2. 保持与现有代码兼容
3. 逐步重构，每步验证编译和功能
4. 最终目标：通过基类接口直接调用 `getRRGivenYTP`

## 编译问题修复历史

### 问题1: const 限定符不匹配
- **文件**: `src/ChemistryModel/FastChemistryModel.C:748`
- **错误**: 函数定义缺少 `const` 限定符
- **解决**: 在函数定义后添加 `const`

```cpp
// 修复前
Foam::scalarField Foam::FastChemistryModel<ThermoType>::getRRGivenYTP(...)

// 修复后
Foam::scalarField Foam::FastChemistryModel<ThermoType>::getRRGivenYTP(...) const
```

## 参考文档
- [化学模型选择机制教程](../docs/ChemistryModelSelection_Tutorial.md)

## 注意事项
1. 编译前确保已正确加载 OpenFOAM 环境
2. 在 `system/controlDict` 中需要加载自定义库：
   ```cpp
   libs
   (
       "libFastChemistryModel.so"
   );
   ```
3. 确保所有函数签名（包括 const 限定符）在声明和定义中完全一致
