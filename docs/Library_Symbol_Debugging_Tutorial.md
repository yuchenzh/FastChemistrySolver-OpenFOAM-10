# C++ 动态库符号调试完整教程

## 📚 目录

1. [什么是 Name Mangling](#什么是-name-mangling)
2. [基础工具介绍](#基础工具介绍)
3. [常用命令速查](#常用命令速查)
4. [实战案例：追踪依赖问题](#实战案例追踪依赖问题)
5. [Name Mangling 编码规则](#name-mangling-编码规则)
6. [常见问题与解决方案](#常见问题与解决方案)

---

## 什么是 Name Mangling

### 概念

**Name Mangling**（名称修饰）是 C++ 编译器将函数名、类名等标识符转换为链接器可以识别的唯一符号的过程。

### 为什么需要 Name Mangling？

C++ 支持以下特性，需要在符号中编码更多信息：

1. **函数重载** - 同名函数但参数不同
2. **命名空间** - 需要包含完整的作用域信息
3. **模板** - 需要包含模板参数类型
4. **类成员函数** - 需要包含类名信息

### 示例对比

```cpp
// C++ 源代码
namespace Foam {
    class FastChemistryModel {
        double solve(double deltaT);
        double solve(const Field<double>& deltaT);
    };
}
```

**编译后的符号**:
```
_ZN4Foam18FastChemistryModel5solveEd                      // solve(double)
_ZN4Foam18FastChemistryModel5solveERKNS_5FieldIdEE        // solve(Field<double> const&)
```

---

## 基础工具介绍

### 1. `nm` - 列出符号表

查看目标文件或库文件中的符号。

**基本语法**:
```bash
nm [选项] <文件>
```

**常用选项**:
- `-D` / `--dynamic`: 只显示动态符号（用于共享库 .so 文件）
- `-C` / `--demangle`: 自动解码 C++ 符号（等同于 `nm | c++filt`）
- `-g`: 只显示外部符号
- `-u`: 只显示未定义的符号

**符号类型**:
| 符号 | 含义 | 说明 |
|------|------|------|
| `T` | Text | 函数代码在此库中已定义（导出） |
| `t` | local text | 本地函数（不导出） |
| `U` | Undefined | 需要从其他库链接的符号 |
| `W` | Weak | 弱符号（可以被覆盖） |
| `D` | Data | 已初始化的数据段 |
| `B` | BSS | 未初始化的数据段 |

### 2. `c++filt` - 解码 C++ 符号

将 mangled 符号转换为可读的 C++ 名称。

**基本语法**:
```bash
c++filt <mangled_name>
echo "<mangled_name>" | c++filt
```

**示例**:
```bash
$ c++filt _ZN4Foam18FastChemistryModel5solveEd
Foam::FastChemistryModel::solve(double)
```

### 3. `ldd` - 列出动态依赖

显示程序或共享库依赖的其他共享库。

**基本语法**:
```bash
ldd <可执行文件或库>
```

**输出格式**:
```
库名 => 实际路径 (加载地址)
```

### 4. `objdump` - 查看目标文件信息

更强大的工具，可以查看反汇编代码、符号表等。

**基本语法**:
```bash
objdump [选项] <文件>
```

**常用选项**:
- `-T`: 显示动态符号表
- `-t`: 显示符号表
- `-d`: 反汇编
- `-C`: 解码 C++ 符号

---

## 常用命令速查

### 📋 基础检查命令

```bash
# 1. 查看库中所有动态符号（未解码）
nm -D libMyLib.so

# 2. 查看库中所有动态符号（解码）
nm -D libMyLib.so | c++filt

# 3. 只显示已定义的函数（T 符号）
nm -D libMyLib.so | grep " T "

# 4. 只显示已定义的函数（解码）
nm -D libMyLib.so | grep " T " | c++filt

# 5. 只显示未定义的符号（可能导致链接错误）
nm -D libMyLib.so | grep " U "

# 6. 只显示未定义的符号（解码）
nm -D libMyLib.so | grep " U " | c++filt
```

### 🔍 搜索特定符号

```bash
# 7. 搜索包含 "solve" 的符号
nm -D libMyLib.so | grep "solve"

# 8. 搜索包含 "solve" 的符号（解码）
nm -D libMyLib.so | grep "solve" | c++filt

# 9. 搜索特定类的符号
nm -D libMyLib.so | grep "FastChemistryModel" | c++filt

# 10. 检查特定符号是否存在
nm -D libMyLib.so | grep "_ZN4Foam18FastChemistryModel5solveEd"
```

### 🔗 检查依赖关系

```bash
# 11. 查看库的所有依赖
ldd libMyLib.so

# 12. 查看依赖并过滤特定库
ldd libMyLib.so | grep "FastChemistry"

# 13. 检查库是否能找到所有依赖（查找 "not found"）
ldd libMyLib.so | grep "not found"
```

### 🛠️ 解码单个符号

```bash
# 14. 解码单个符号
echo "_ZN4Foam18FastChemistryModel5solveEd" | c++filt

# 15. 直接用 c++filt
c++filt _ZN4Foam18FastChemistryModel5solveEd

# 16. 批量解码多个符号
cat symbols.txt | c++filt
```

### 📊 高级分析

```bash
# 17. 统计库中函数数量
nm -D libMyLib.so | grep " T " | wc -l

# 18. 列出所有类
nm -D libMyLib.so | c++filt | grep "::" | awk -F'::' '{print $1}' | sort -u

# 19. 查找未定义的 OpenFOAM 符号
nm -D libMyLib.so | grep " U " | grep "Foam" | c++filt

# 20. 比较两个版本的符号差异
diff <(nm -D libOld.so | sort) <(nm -D libNew.so | sort)
```

### 🐛 调试命令

```bash
# 21. 查看编译时链接了哪些库（从 Make/options）
grep -E "^LIB_LIBS" Make/options

# 22. 检查 OpenFOAM 环境变量
echo $FOAM_USER_LIBBIN
ls -la $FOAM_USER_LIBBIN/*.so

# 23. 查看库的导出符号数量统计
nm -D libMyLib.so | awk '{print $2}' | sort | uniq -c

# 24. 使用 objdump 查看动态符号
objdump -T libMyLib.so | grep "solve"

# 25. 检查库的架构和格式
file libMyLib.so
```

---

## 实战案例：追踪依赖问题

### 问题场景

运行 `reactingFoam` 时出现错误：
```
dlopen error : libFastChemistryModel.so: undefined symbol: _ZN4Foam18FastChemistryModel5solveEd
--> FOAM Warning : could not load "libCCM.so"
```

### 调试步骤

#### 第 1 步：解码错误信息

```bash
$ c++filt _ZN4Foam18FastChemistryModel5solveEd
Foam::FastChemistryModel::solve(double)
```

**结论**: 缺少 `FastChemistryModel::solve(double)` 函数的实现。

#### 第 2 步：检查库中的符号

```bash
$ nm -D $FOAM_USER_LIBBIN/libFastChemistryModel.so | grep "FastChemistryModel.*solve"
                 U _ZN4Foam18FastChemistryModel5solveEd
```

**发现**: 符号类型是 `U` (Undefined)，说明库声明了这个函数但没有提供实现！

#### 第 3 步：查找函数声明

```bash
$ grep -rn "solve.*deltaT" src/ChemistryModel/FastChemistryModel.H
163:        scalar solve(const DeltaTType& deltaT){return 0;}
```

**发现**: 头文件中有声明。

#### 第 4 步：查找函数实现

```bash
$ find src -name "*solve*.H"
src/ChemistryModel/FastChemistryModel_transientSolve.H
src/ChemistryModel/FastChemistryModel_localEulerSolve.H
```

**发现**: 实现在单独的 `.H` 文件中。

#### 第 5 步：检查是否包含实现

```bash
$ grep -n "#include.*solve" src/ChemistryModel/FastChemistryModel.C
# 没有输出
```

**根本原因**: `FastChemistryModel.C` 没有 `#include` 这些实现文件！

#### 第 6 步：修复并验证

添加 include：
```cpp
// FastChemistryModel.C 末尾
#include "FastChemistryModel_transientSolve.H"
#include "FastChemistryModel_localEulerSolve.H"
```

重新编译后检查：
```bash
$ wmake -j
$ nm -D $FOAM_USER_LIBBIN/libFastChemistryModel.so | grep "FastChemistryModel.*solve"
000000000002a890 T _ZN4Foam18FastChemistryModel5solveEd
000000000002a880 T _ZN4Foam18FastChemistryModel5solveERKNS_5FieldIdEE
```

**成功**: 符号类型从 `U` 变成了 `T`！

#### 第 7 步：验证依赖加载

```bash
$ ldd $FOAM_USER_LIBBIN/libCCM.so | grep Fast
libFastChemistryModel.so => /path/to/libFastChemistryModel.so (0x...)
```

**成功**: CCM 库能找到 FastChemistryModel 库。

#### 第 8 步：测试运行

```bash
$ cd Sandia
$ reactingFoam
# 不再出现 dlopen error!
```

---

## Name Mangling 编码规则

### 基本结构

```
_Z                        → C++ mangling 魔法前缀
N...E                     → 嵌套名称（类/命名空间）
<长度><名称>               → 长度编码的标识符
```

### 类型编码表

#### 基本类型

| Mangled | C++ 类型 | 说明 |
|---------|---------|------|
| `v` | `void` | 无返回值 |
| `b` | `bool` | 布尔值 |
| `c` | `char` | 字符 |
| `i` | `int` | 整数 |
| `j` | `unsigned int` | 无符号整数 |
| `l` | `long` | 长整数 |
| `m` | `unsigned long` | 无符号长整数 |
| `x` | `long long` | 64位整数 |
| `y` | `unsigned long long` | 无符号64位整数 |
| `f` | `float` | 单精度浮点 |
| `d` | `double` | 双精度浮点 |

#### 修饰符

| Mangled | C++ 含义 | 说明 |
|---------|---------|------|
| `P` | `*` | 指针 (Pointer) |
| `R` | `&` | 引用 (Reference) |
| `K` | `const` | 常量修饰 |
| `V` | `volatile` | 易变修饰 |
| `r` | `restrict` | 限制修饰 |

#### 特殊结构

| Mangled | 含义 | 说明 |
|---------|------|------|
| `N...E` | 嵌套名称 | 命名空间::类::函数 |
| `I...E` | 模板参数 | 模板类型列表 |
| `S_` | 替换 | 重用之前的类型/命名空间 |
| `St` | `std::` | 标准库命名空间 |

### 解析示例

#### 示例 1：简单函数

```
_ZN4Foam18FastChemistryModel5solveEd

分解：
_Z                          → mangling 前缀
N                           → 开始嵌套名称
  4Foam                     → 命名空间 "Foam" (长度4)
  18FastChemistryModel      → 类名 "FastChemistryModel" (长度18)
  5solve                    → 函数名 "solve" (长度5)
E                           → 结束嵌套名称
d                           → 参数: double

结果：Foam::FastChemistryModel::solve(double)
```

#### 示例 2：带引用参数

```
_ZN4Foam18FastChemistryModel5solveERKNS_5FieldIdEE

分解：
_Z                          → mangling 前缀
N                           → 开始嵌套名称
  4Foam                     → 命名空间 "Foam"
  18FastChemistryModel      → 类名
  5solve                    → 函数名
E                           → 结束嵌套名称

RK                          → Reference to Const (const &)
NS_                         → Nested, 重用命名空间 (S_ = Foam)
  5Field                    → 类名 "Field"
  I                         → 开始模板参数
    d                       → double
  E                         → 结束模板参数
E                           → 结束整个类型

结果：Foam::FastChemistryModel::solve(Foam::Field<double> const&)
```

#### 示例 3：指针参数

```
_ZN8LUsolver11printMatrixEPdii

分解：
_Z                          → mangling 前缀
N                           → 开始嵌套名称
  8LUsolver                 → 类名 "LUsolver" (长度8)
  11printMatrix             → 函数名 "printMatrix" (长度11)
E                           → 结束嵌套名称

Pd                          → Pointer to double (double*)
i                           → int
i                           → int

结果：LUsolver::printMatrix(double*, int, int)
```

#### 示例 4：模板类

```
_ZNK4Foam10OptRodas34INS_18FastChemistryModelEE5solveEidRdS3_

分解：
_Z                          → mangling 前缀
NK                          → Nested + Const (const 成员函数)
  4Foam                     → 命名空间 "Foam"
  10OptRodas34              → 类名 "OptRodas34"
  I                         → 开始模板参数
    NS_18FastChemistryModelE → Foam::FastChemistryModel
  E                         → 结束模板参数
  5solve                    → 函数名 "solve"
E                           → 结束嵌套名称

i                           → int
d                           → double
Rd                          → double& (引用)
S3_                         → 替换第3个类型 (double&)

结果：Foam::OptRodas34<Foam::FastChemistryModel>::solve(int, double, double&, double&) const
```

### 手动验证长度

```bash
$ echo -n "Foam" | wc -c
4

$ echo -n "FastChemistryModel" | wc -c
18

$ echo -n "solve" | wc -c
5
```

---

## 常见问题与解决方案

### 问题 1：undefined symbol 错误

**错误信息**:
```
dlopen error: libMyLib.so: undefined symbol: _ZN...
```

**可能原因**:
1. 函数声明了但没有定义
2. 链接时缺少依赖库
3. 模板函数没有实例化

**调试步骤**:
```bash
# 1. 解码符号
echo "符号名" | c++filt

# 2. 检查库中是否有这个符号
nm -D libMyLib.so | grep "符号名"

# 3. 检查符号类型
# U = 未定义（问题）
# T = 已定义（正常）

# 4. 检查依赖库
ldd libMyLib.so

# 5. 查找函数实现
grep -r "函数名" src/
```

### 问题 2：符号冲突

**错误信息**:
```
multiple definition of `Foo::bar()'
```

**可能原因**:
1. 同一个函数在多个 .C 文件中定义
2. 头文件中的非 inline 函数被多次包含
3. 模板特化冲突

**解决方案**:
- 使用 `inline` 关键字
- 使用匿名命名空间
- 检查重复包含

### 问题 3：ABI 不兼容

**错误信息**:
```
version `GLIBCXX_3.4.26' not found
```

**可能原因**:
编译库的编译器版本与运行环境的编译器版本不匹配。

**检查方法**:
```bash
# 查看库需要的 GLIBCXX 版本
strings libMyLib.so | grep GLIBCXX

# 查看系统提供的版本
strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX
```

### 问题 4：找不到共享库

**错误信息**:
```
error while loading shared libraries: libMyLib.so: cannot open shared object file
```

**解决方案**:
```bash
# 1. 检查库路径
echo $LD_LIBRARY_PATH

# 2. 添加库路径
export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH

# 3. 或在 controlDict 中加载
libs ("libMyLib.so");
```

---

## 实用脚本

### 符号检查脚本

创建 `check_symbols.sh`:
```bash
#!/bin/bash
# 用法: ./check_symbols.sh libMyLib.so [搜索关键词]

LIB=$1
KEYWORD=${2:-""}

if [ ! -f "$LIB" ]; then
    echo "错误: 库文件 $LIB 不存在"
    exit 1
fi

echo "========================================="
echo "库文件: $LIB"
echo "========================================="

if [ -n "$KEYWORD" ]; then
    echo "搜索关键词: $KEYWORD"
    echo ""
    echo "已定义的函数 (T):"
    nm -D "$LIB" | grep " T " | grep "$KEYWORD" | c++filt
    echo ""
    echo "未定义的符号 (U):"
    nm -D "$LIB" | grep " U " | grep "$KEYWORD" | c++filt
else
    echo "符号统计:"
    echo "  已定义函数 (T): $(nm -D "$LIB" | grep " T " | wc -l)"
    echo "  未定义符号 (U): $(nm -D "$LIB" | grep " U " | wc -l)"
    echo ""
    echo "前 10 个导出函数:"
    nm -D "$LIB" | grep " T " | c++filt | head -10
fi
```

使用：
```bash
chmod +x check_symbols.sh
./check_symbols.sh libFastChemistryModel.so solve
```

### 依赖检查脚本

创建 `check_deps.sh`:
```bash
#!/bin/bash
# 用法: ./check_deps.sh libMyLib.so

LIB=$1

if [ ! -f "$LIB" ]; then
    echo "错误: 库文件 $LIB 不存在"
    exit 1
fi

echo "========================================="
echo "库依赖检查: $LIB"
echo "========================================="

echo ""
echo "直接依赖:"
ldd "$LIB" | grep "=>"

echo ""
echo "缺失的依赖:"
MISSING=$(ldd "$LIB" | grep "not found")
if [ -z "$MISSING" ]; then
    echo "  无"
else
    echo "$MISSING"
fi

echo ""
echo "OpenFOAM 库依赖:"
ldd "$LIB" | grep OpenFOAM
```

---

## 快速参考卡片

```
┌─────────────────────────────────────────────────────┐
│  C++ 动态库调试速查卡                                 │
├─────────────────────────────────────────────────────┤
│ 查看符号:                                            │
│   nm -D lib.so                    原始符号           │
│   nm -D lib.so | c++filt          解码符号           │
│                                                      │
│ 搜索符号:                                            │
│   nm -D lib.so | grep "solve" | c++filt              │
│                                                      │
│ 符号类型:                                            │
│   T = 已定义 (正常)    U = 未定义 (可能有问题)        │
│   W = 弱符号           t = 本地函数                   │
│                                                      │
│ 解码单个符号:                                         │
│   echo "_ZN..." | c++filt                            │
│   c++filt _ZN...                                     │
│                                                      │
│ 检查依赖:                                            │
│   ldd lib.so                      所有依赖           │
│   ldd lib.so | grep "not found"   缺失依赖           │
│                                                      │
│ 调试流程:                                            │
│   1. 解码错误中的符号                                 │
│   2. 检查库中符号类型 (T/U)                           │
│   3. 查找函数声明和定义                               │
│   4. 检查编译配置 (Make/files, Make/options)         │
│   5. 修复并重新编译                                   │
│   6. 验证符号从 U 变成 T                              │
└─────────────────────────────────────────────────────┘
```

---

## 进阶话题

### 1. 使用 `readelf` 查看 ELF 信息

```bash
# 查看动态符号表
readelf -s --dyn-syms libMyLib.so

# 查看依赖库
readelf -d libMyLib.so | grep NEEDED
```

### 2. 使用 `objdump` 反汇编

```bash
# 查看特定函数的汇编代码
objdump -d libMyLib.so | less

# 查找特定函数
objdump -d libMyLib.so | grep -A 20 "solve"
```

### 3. 使用 `strings` 查看字符串

```bash
# 查看库中的所有字符串
strings libMyLib.so | less

# 查找特定字符串
strings libMyLib.so | grep "error"
```

---

## 参考资料

### 在线工具

- **C++ Name Demangler**: https://demangler.com/
- **GCC Name Mangling**: https://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html

### 相关文档

- `man nm` - nm 命令手册
- `man c++filt` - c++filt 命令手册
- `man ldd` - ldd 命令手册
- `man objdump` - objdump 命令手册

### OpenFOAM 特定

- OpenFOAM Wiki: https://openfoamwiki.net/
- OpenFOAM Source: https://github.com/OpenFOAM/OpenFOAM-dev

---

## 总结

掌握这些工具和技巧，你就能：

✅ 快速解码神秘的符号错误
✅ 追踪库依赖问题的根本原因
✅ 理解 C++ 编译和链接过程
✅ 独立调试动态库加载问题

记住核心流程：
1. **解码** → `c++filt`
2. **检查** → `nm -D`
3. **分析** → `ldd`
4. **修复** → 添加实现/依赖
5. **验证** → 重新编译检查

---

**文档版本**: 1.0
**创建日期**: 2025-11-05
**作者**: Claude Code Tutorial
**基于案例**: FastChemistrySolver-OpenFOAM-10 依赖问题调试
