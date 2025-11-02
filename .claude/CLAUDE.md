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

### 目标
将 `FastChemistryModel` 的 `getRRGivenYTP` 功能包装为一个通用接口，使得其他化学模型也能轻松使用该功能，而无需显式修改 `basicChemistryModel` 基类。

### 设计约束
1. 不直接修改 `basicChemistryModel` 基类添加 `getRRGivenYTP` 虚函数
2. 提供一个灵活的包装器或工具类
3. 保持与现有代码的兼容性
4. 易于在 `testFoam` 等求解器中调用

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
