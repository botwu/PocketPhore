# PhoreGen原始数据集详细分析

## 1. 数据集概览

PhoreGen使用了三个数据集，**均只包含药效团信息，无蛋白质口袋数据**：

| 数据集 | 来源 | 数量 | 内容 | 下载地址 |
|--------|------|------|------|---------|
| **LigPhore** | ZINC 300万分子 | ~训练集数十万 | 配体 + 药效团 | Zenodo: 15518867 |
| **CpxPhore** | PDBbind复合物 | ~数千 | 配体 + 药效团(从复合物提取) | Zenodo: 15518867 |
| **DockPhore** | CrossDocked对接 | ~数千 | 配体 + 药效团(从对接pose提取) | Zenodo: 15518867 |

## 2. 数据格式

### .phore文件格式
```
ZINC_ID_conformer_featuretype
HA	1.000	1.200	1.000	0.210	-1.197	-0.159	1	1.915	-1.543	-1.663	6	1.000
AR	0.700	1.000	1.000	-3.627	-0.909	-0.238	1	-4.053	2.024	-0.703	0	1.000
HD	1.000	1.200	1.000	1.363	2.438	-1.625	1	-0.143	4.036	-2.310	5	1.000
...
```

**字段含义**：
- Column 1: 特征类型 (HA/HD/AR/HY/PO/NE/CV/XB/EX)
- Column 2-4: alpha, weight, factor (权重参数)
- Column 5-7: x, y, z (3D坐标)
- Column 8: has_norm (是否有方向向量)
- Column 9-11: norm_x, norm_y, norm_z (方向向量)
- Column 12: label
- Column 13: anchor_weight

### 药效团特征类型
```python
PHORETYPES = ['MB', 'HD', 'AR', 'PO', 'HA', 'HY', 'NE', 'CV', 'CR', 'XB', 'EX']

# MB: Metal Binding
# HD: Hydrogen Bond Donor
# AR: Aromatic Ring
# PO: Positive Ionizable
# HA: Hydrogen Bond Acceptor
# HY: Hydrophobic
# NE: Negative Ionizable
# CV: Covalent (分为CV1-CV4四类)
# CR: Chelating Ring (训练时跳过)
# XB: Halogen Bond
# EX: Exclusion Volume (排斥体积)
```

## 3. 关键发现：PhoreGen缺少的信息

### ❌ 没有的信息
1. **蛋白质结构** - 无PDB文件
2. **口袋残基** - 无残基序列和坐标
3. **蛋白-配体相互作用** - 无氢键、疏水等信息
4. **结合亲和力** - 虽然用了PDBbind，但未直接用于训练

### ✅ 有的信息
1. **药效团模型** - 位置、类型、方向
2. **配体结构** - SMILES、3D坐标、原子类型、键类型
3. **药效团方向信息** - Direction Matching用到的norm向量

## 4. 数据生成流程

```
ZINC配体分子
    ↓
使用LigandScout/RDKit
    ↓
提取药效团特征 (HA, HD, AR, HY等)
    ↓
保存为.phore文件
    ↓
配对配体分子数据
    ↓
LigPhore数据集
```

对于CpxPhore/DockPhore：
```
PDBbind蛋白-配体复合物
    ↓
提取配体分子
    ↓
使用LigandScout (structure-based模式)
    ↓
考虑蛋白环境生成药效团
    ↓
保存.phore文件 (但不保存蛋白)
    ↓
CpxPhore/DockPhore数据集
```

## 5. 训练策略

### 两阶段训练
```python
# 阶段1: LigPhore预训练
python train.py --config ./configs/train_lig-phore.yml
# 训练160 epochs，学习基础药效团→配体生成

# 阶段2: CpxPhore + DockPhore微调
python train.py --config ./configs/train_dock-cpx-phore.yml
# 微调学习更真实的蛋白环境中的药效团
```

### 数据增强
```yaml
add_lig_noise: True
lig_noise_std: 0.1          # 配体坐标加噪声
add_phore_noise: True
phore_noise_std: 0.1        # 药效团坐标加噪声
phore_norm_angle: 5.0       # 药效团方向扰动5度
```

## 6. 关键配置

```yaml
# 从 train_lig-phore.yml
dataset:
  zinc_train_filelist: ./data/index/dense/train_filelist.pkl
  zinc_valid_filelist: ./data/index/dense/valid_filelist.pkl
  zinc_test_filelist: ./data/index/dense/test_filelist.pkl
  save_path: ./data/ZINC_300W_dense_pkl
  data_name: zinc_300
  max_atom: 78
  center: phore  # 以药效团为中心对齐

model:
  phore_feat_dim: 18  # 药效团特征维度
  # 对于zinc_300和pdbbind，会自动+2变为18
  # 包含: one-hot(11) + alpha(1) + has_norm(2) + exclusion(2) = 16 → 18
```

## 7. 为什么PhoreGen不需要蛋白质信息？

PhoreGen的设计哲学是：
- **药效团已经隐含了蛋白质约束**
  - CpxPhore的药效团是从蛋白-配体复合物中提取的
  - 药效团位置和类型已经编码了关键的相互作用信息

- **关注生成多样性而非特定结合能力**
  - 目标是生成满足药效团约束的分子
  - 不强调与特定蛋白口袋的精确匹配

- **简化训练复杂度**
  - 避免处理大量蛋白质结构数据
  - 降低数据收集和预处理成本

## 8. PhoreGen的局限性（也是我们改进的动机）

| 局限性 | 后果 | 改进方向 |
|--------|------|---------|
| 无蛋白质几何约束 | 可能生成与口袋冲突的分子 | 加入口袋残基信息 |
| 无物理化学约束 | 氢键/疏水等相互作用不精确 | 物理约束模块 |
| 药效团信息静态 | 无法动态适应不同口袋大小 | 口袋感知的大小预测 |
| 单一尺度引导 | 缺乏全局/局部层次化优化 | 双尺度注意力机制 |

这正是我们提出PocketPhore改进方案的原因！
