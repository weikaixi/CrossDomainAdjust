# CrossDomainAdjust：跨域特征转换与部分校正

**版本 2.1.0。** [English](README.md) · [中文完整使用手册](inst/doc/USER_GUIDE_zh.md) ·
[英文完整使用手册](inst/doc/USER_GUIDE.md) · [15 个可运行案例](inst/examples)

这个 R 包在训练集学习转换规则，将特征集合、配对筛选、标准化参数、域中心和投影方向保存下来，
然后应用于外部队列或单个新样本。域可以表示研究中心、技术平台，也可以表示 LUAD/LUSC、LGG/GBM
等生物学类别。生物学域差异不等同于实验批次效应，因此包提供部分校正，而不要求全部消除。

| 方法 | 参数 | 实际处理 |
|---|---|---|
| 排序 | `method="rank"` | 样本内分数排序 → 训练集冻结 Z-score |
| 配对 | `method="pair"` | 两基因比较得到 0/1 → 训练集频率 20%–80% 筛选 → 冻结 Z-score |
| 域中心平移 | `method="center"` | 原始输入表达减去 λ 倍的域中心偏移，不额外加 Z-score |
| 域投影 | `method="project"` | 原始输入表达减去 λ 倍的域方向分量，不额外加 Z-score |

“原始输入表达”指已经完成适当平台预处理的表达矩阵，不是未经处理的 RNA-seq counts。
λ 默认候选值为 0、0.2、0.4、0.6、0.8、1；排序和配对本身不使用 λ。

## 安装

需要 R >= 4.1.0；本包不含编译代码，不需要为本包安装 Rtools。

```r
install.packages("remotes")
remotes::install_github("weikaixi/CrossDomainAdjust", upgrade = "never")
library(CrossDomainAdjust)
packageVersion("CrossDomainAdjust")
```

## 最小完整案例

```r
d <- simulate_domain_data(p = 8)
fit <- fit_domain_method(d$train$x, d$train$domain, "project", lambda = 0.4)
train_new <- transform_cross_domain(fit, d$train$x)
test_new <- transform_cross_domain(fit, d$test$x)
single_new <- transform_cross_domain(fit, d$test$x[1, ])
lambda_diagnostics(fit, d$train$x, d$train$domain)$summary
```

中心平移在变换时还需提供域标签：`transform_cross_domain(fit, x, domain)`。
配对可以加 `frequency=c(0.2,0.8), max_pairs=128`，控制保留的特征数量。
后续接 Cox、逻辑回归、线性回归等模型；本包负责跨域预处理，不把某个预测算法固定在包里。

## 15 个案例

1. 四种方法的统一入口。
2. 排序的单调变换不变性及冻结标准化。
3. 配对频率筛选、显式候选配对和特征数量上限。
4. 已知域的中心平移及未知域报错。
5. 未知来源的单个样本投影。
6. λ 网格、固定 PCA 坐标和一致坐标范围。
7. 回归任务中的折内调参和独立测试。
8. 二分类与 Brier 分数。
9. 生存分析、Ridge Cox 与 pooled C-index。
10. 多域、稳健中心与投影维数。
11. 多种子敏感性分析，保留完整网格。
12. RDS 保存、加载、乱序基因和缺失基因校验。
13. 高级组合接口：中心平移后投影。
14. 真实肺癌：TCGA LUAD/LUSC → GSE3141。
15. 真实胶质瘤：TCGA LGG/GBM → GSE43378。

```r
source(system.file("examples", "14_real_lung_survival.R", package="CrossDomainAdjust"))
source(system.file("examples", "15_real_glioma_survival.R", package="CrossDomainAdjust"))
source(system.file("examples", "run_all.R", package="CrossDomainAdjust"))
```

前 13 个使用模拟教学数据；后 2 个使用随包提供的公开来源、经过处理的 60 基因数据。
真实数据案例为了演示接口采用固定 λ 与惩罚强度，不能把其输出误认为论文冻结模型的原始结果。

λ=1 可以使拟合数据的域中心重合，但不代表协方差、所有分布差异或外部批次效应都消失。
对预后有帮助的域相关生物学信号也可能被消除。λ 应由训练阶段验证确定，不强制选在 0.2–0.8。
完整公式、参数解释、真实数据准备、调参细节、常见错误与部署步骤见[中文手册](inst/doc/USER_GUIDE_zh.md)。
