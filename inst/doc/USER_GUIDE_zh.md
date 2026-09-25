# CrossDomainAdjust 2.1.0 中文完整使用手册

## 1. 这个包解决什么问题

本包把“训练时学到的预处理规则”冻结下来，再应用到另一批样本，避免新队列到来后重新定义特征、
重新筛选配对或重新估计标准化参数。核心有四种方法：排序、配对、域中心平移、域投影。
前两种主要改变特征表示；后两种直接对输入表达做部分域校正。

域不一定是实验批次。肺癌中的 LUAD/LUSC、胶质瘤中的 LGG/GBM 是生物学类别；中心、平台、
实验室则可以是技术分组。一个域方向可能同时包含组织学、生物学和技术信息。
因此本包不假设所有域差异都应被消除，也不承诺某种方法在所有队列中最好。

本包不自动完成测序 counts 归一化、探针注释、基因 ID 映射、缺失值填补、预后基因选择或临床终点清洗。
这些步骤需在进入本包前确定；涉及训练参数或结局选择时，必须限制在相应训练数据中。

## 2. 安装和查看帮助

```r
install.packages("remotes")
remotes::install_github("weikaixi/CrossDomainAdjust", upgrade="never")
library(CrossDomainAdjust)
packageVersion("CrossDomainAdjust")
help(package="CrossDomainAdjust")
?fit_domain_method
?tune_cross_domain
```

R 版本至少 4.1.0。核心只依赖 R 自带的 stats、utils，不需要为本包编译 C/C++。
生存案例使用 survival，通常随 R 一起安装；缺少时运行 `install.packages("survival")`。
本地源码安装：`install.packages("CrossDomainAdjust_2.1.0.tar.gz", repos=NULL, type="source")`。

安装后，`system.file("examples", package="CrossDomainAdjust")` 是案例目录，
`system.file("doc", package="CrossDomainAdjust")` 是中英文手册目录。
长期复现时除了包版本，还应记录 Git 提交号或源码压缩包校验值，避免默认分支更新后难以追溯。

## 3. 输入矩阵怎么整理

**行是样本，列是基因。** 列名唯一且不能为空，数值必须有限，不能有 NA、NaN、Inf。
一个新样本可以使用带基因名称的数值向量。数值型 data.frame 也可以，但不能把文字样本 ID 放在数值列中。
当前不直接接受稀疏矩阵和延迟矩阵。

```r
# 常见基因×样本表的读取方式，取消注释并替换路径后使用：
# expr <- read.delim("expression.tsv", row.names=1, check.names=FALSE)
# x <- t(as.matrix(expr))
# cl <- read.csv("clinical.csv")
# stopifnot(!anyDuplicated(cl$sample_id))
# cl <- cl[match(rownames(x), cl$sample_id), ]
# stopifnot(!anyNA(cl$sample_id), identical(rownames(x), cl$sample_id))
# domain <- cl$histology
```

训练与预测必须提供同一套已冻结基因，列顺序可以不同，程序会按名字对齐。额外的有限数值列会被忽略，
但缺少任何拟合时的基因都会报错。即使配对后只使用部分基因，当前实现仍要求完整的拟合基因集合。
不能把外部缺失的基因直接当成 0，也不能不说明就改为另一套基因重新排名。

如果原始数据有多个探针对应一个基因，要事先规定合并方法；不能根据外部生存结果选择“最好探针”。
需要填补缺失值时，应由训练数据学习填补规则，而不是让测试集改变训练规则。
原始表达指已经适当处理的表达尺度，不是未经处理的 RNA-seq 原始计数。

## 4. 四种方法的统一入口

```r
d <- simulate_domain_data(p=8, seed=42)
tr <- d$train
te <- d$test

rank_fit <- fit_domain_method(tr$x, tr$domain, "rank")
pair_fit <- fit_domain_method(tr$x, tr$domain, "pair", frequency=c(.2,.8), max_pairs=128)
center_fit <- fit_domain_method(tr$x, tr$domain, "center", lambda=.4)
project_fit <- fit_domain_method(tr$x, tr$domain, "project", lambda=.4)

rank_new <- transform_cross_domain(rank_fit, te$x)
pair_new <- transform_cross_domain(pair_fit, te$x)
center_new <- transform_cross_domain(center_fit, te$x, te$domain)
project_new <- transform_cross_domain(project_fit, te$x)
```

排序和配对预设都包含训练集冻结 Z-score；中心平移和域投影使用原始输入尺度，不额外加 Z-score。
λ 只控制平移和投影强度，排序、配对预设不使用 λ，传入也会被忽略。

## 5. 排序为什么可能有用，以及排序后是否标准化

在每个样本内部，对 p 个基因做排序，转换成 `(rank-1)/(p-1)`，平秩采用平均名次。
这个步骤使分数排名落在 0–1 之间。它不是在所有样本之间对一个基因排名。
然后，rank 预设再用训练集每个排名特征的均值和样本标准差做 Z-score。
两者作用不同：前者处理样本内相对表达，后者让进入对尺度敏感模型的特征更可比。

一个样本所有基因同时加相同常数、乘正数，或受到相同严格递增变换时，其基因顺序不变。
因此排名对某些整体尺度和位置变化不敏感。但每个基因受不同的平台偏移、存在大量平秩或缺失基因时，
这个性质不能保证排序保持不变。排序好用有数学依据，但不能推出对任何批次都有效。

```r
z1 <- transform_cross_domain(rank_fit, te$x)
z2 <- transform_cross_domain(rank_fit, exp(te$x/10))
stopifnot(isTRUE(all.equal(z1,z2)))
range(encode_features(rank_fit$encoder,te$x)) # 原始分数排名
range(z1)                                  # Z-score 后可以超出 0–1
```

在外部数据上不能另外运行 `scale(test)` 重新估计均值和标准差。冻结后同一新样本单独输入，
应与其在整批样本中转换的结果一致。训练常数特征的缩放分母设为 1；不会自动删除这些特征。
标准化不是对所有模型都必然更好。高级接口允许 rank-only，但论文四方法预设采用冻结 Z-score。

## 6. 配对特征如何筛选，如何避免数量爆炸

对于基因 i、j，计算 `I(x_i>x_j)`。大于为 1，小于或等于为 0。
每个无序基因对只比较一次，方向由训练矩阵列顺序确定。显式输入反向配对也会统一到该顺序；
重复或正反方向重复的候选会报错。

在训练集统计 1 的比例，默认保留 **0.2≤频率≤0.8**，端点包含在内。
这是训练集总体频率筛选，不是每一个域分别都满足频率要求；大样本域会对总体频率影响更大。
筛选后再冻结 Z-score，最终输入模型的值不再是严格的 0/1。

p 个基因有 p(p−1)/2 个候选对。60 个基因只有 1770 个候选；2 万个基因则接近 2 亿个。
频率筛选减少输出数量，但仍要计算候选，所以不能仅依靠最后筛选解决无限制枚举的问题。

| 参数 | 含义 |
|---|---|
| `frequency=c(.2,.8)` | 保留训练集 1 的出现比例范围 |
| `max_pairs=128` | 频率筛选后最多保留多少个特征 |
| `max_candidates=1000000` | 枚举前候选数量上限，超过直接报错 |
| `chunk_size=1000` | 分块计算候选比较，降低中间矩阵内存 |
| `pairs=` | 显式指定两列基因名，避免枚举所有配对 |

达到 `max_pairs` 上限时，优先保留频率接近 0.5 的配对，相同时按候选顺序。
这不是按预后 P 值筛选，也不等同于去除共线性；冗余配对仍可能存在。
如果没有候选通过筛选，程序报错，不会偷偷放宽阈值。

```r
candidates <- rbind(c("gene1","gene2"),c("gene3","gene4"))
fit <- fit_domain_method(tr$x,tr$domain,"pair",pairs=candidates,frequency=c(0,1))
fit$encoder$pairs
fit$encoder$frequencies
```

面对高维矩阵，应首先在开发数据中确定合理的基因面板或候选配对。
`max_pairs` 只限制最后输出，不能替代 `max_candidates`；单纯调大后者可能带来巨大计算成本。

## 7. 中心平移和域投影的公式

设样本为行向量 x，域中心为 μ_d，参照中心为 a，偏移为 δ_d=μ_d−a。

中心平移：`x*=x−λδ_d`。默认使用均值时，拟合数据的域均值变为 `a+(1−λ)δ_d`。
因此 λ 从 0 到 1 时，域中心逐步靠近，并在 1 时重合。
同一域内所有样本减去相同常量，所以域内协方差不会改变。
转换新样本必须知道它属于哪个已拟合域，不能悄悄用外部数据估计新中心。

域投影：把所有 δ_d 堆叠，用奇异值分解求正交方向矩阵 B，再计算 `x*=x−λ(x−a)BB'`。
它把沿着 B 的分量缩小为原来的 1−λ，其余正交方向保留。
程序采用低秩矩阵乘法，不显式创建庞大的基因数×基因数投影矩阵。

投影不需要新样本的域标签，也不需要同时输入一个完整外部队列。
但新平台如果产生了不在已学习方向内的偏移，这个偏移不会自动消失。
与平移不同，投影会改变域内沿 B 方向的变异；有用生物学和技术差异可能一起被削弱。

```r
one <- transform_cross_domain(project_fit,te$x[1,])
batch <- transform_cross_domain(project_fit,te$x)
stopifnot(isTRUE(all.equal(as.numeric(one),as.numeric(batch[1,]))))
```

## 8. 域怎么定义，参照中心怎么选

肺癌真实案例的域标签是 LUAD/LUSC，胶质瘤是 LGG/GBM，不能说它们就是技术批次。
来自 GEO 的 LUAD 可使用 TCGA 学到的 LUAD 中心，因为类别相同；这不代表 GEO 平台偏移已经被单独估计。
如果研究目标是去除技术批次，则应提供技术批次标签；完全未知的批次不能直接使用非零强度中心平移。

默认 `anchor="domain_mean"`，给每个域的中心相同权重。
`anchor="sample_mean"` 根据各域样本量加权；使用均值中心时，它等于所有训练样本的均值。
也可以提供命名的自定义参照向量，坐标尺度必须与送入校正器的特征一致。

默认 `center="mean"`，还可以用 `"median"` 或 `"trimmed"`；截尾均值的 `trim=.1` 表示每端截去 10%。
稳健中心可降低极端值影响，但不保证 λ=1 时算术均值对齐。
在稳健中心下，sample_mean 是稳健中心的样本量加权平均，不再等同原始样本算术均值。

K 个域、内置参照时，投影秩最多 K−1。默认保留所有非零方向；`n_components=1` 等参数可以保留前几个方向。
`n_components=0` 不去除任何投影方向。只保留部分方向时，λ=1 不保证所有中心都重合。
自定义参照位于域中心仿射空间之外时，秩可以达到 K。

## 9. λ 怎么选择

推荐用开发数据内的交叉验证，默认网格 `c(0,.2,.4,.6,.8,1)`。
不要因为想要部分校正的结果，就先删除 0 和 1；端点也应该公平参与选择。
内部最优值可能是 .4，也可能是 0 或 1。看起来两个域完全重合，不等于预测最好。

```r
fit_model <- function(x,y) {
  z <- cbind(1,x)
  solve(crossprod(z)+diag(c(0,rep(1,ncol(x)))),crossprod(z,y))
}
predict_model <- function(model,x) as.numeric(cbind(1,x)%*%model)
score <- function(y,prediction,domain) -mean((y-prediction)^2)
set.seed(7)
folds <- sample(rep(1:3,length.out=nrow(tr$x)))
tuned <- tune_cross_domain(tr$x,tr$domain,tr$regression,folds,
  fit_model,predict_model,score,correction="project")
tuned$best_lambda
tuned$scores
tuned$fold_scores
```

每个训练折单独学习配对筛选、标准化、域中心和投影方向，然后应用于其验证折。
每个 λ 都重新拟合下游模型。函数认为 **分数越大越好**，因此误差指标需要取负。
各折得分等权平均；完全相同的分数选更小的 λ。模型出错或评分非有限值时会停止，不会删除失败折。

函数返回的是 `fit` 预处理对象，不含最终预测模型。选好以后，需要对整个开发集变换，再重新拟合预测模型。
内部选择得分不能作为无偏最终表现，要使用外部未参与选择的数据或嵌套外层验证。
同一患者的多个样本必须放在同一折中。

中心平移要求每个验证折里的域在对应训练折中已出现。留一域验证可使用投影，但不能通过验证数据偷偷估计中心。
如果输入基因是根据全部数据的预后筛选的，这个函数不会自动撤销前面的泄漏；监督式基因筛选也必须进入训练折。
本函数不自动联合选择 λ 和 Cox 惩罚参数。需要联合调参时应明确构建嵌套或二维网格，并保存完整记录。

## 10. 分类、回归、生存模型怎么接入

包只负责特征处理，输出数值矩阵，可以接线性回归、逻辑回归、Cox 或用户自己的模型。
训练回调接收 `(x,y)`，预测回调接收 `(model,x)`，评分回调接收 `(y,prediction,domain)`。
y 可以是向量、矩阵、data.frame 或 survival::Surv，结构化结局会按行切分而保留维度。

回归可以用负 MSE；分类可以用负 Brier 分数或明确定义的 AUC。
分类阈值必须事先确定或在训练阶段选择，不能在外部队列中选择最好的阈值再报告精度。
第 08 个案例使用固定 0.5 阈值，同时报告 Brier 分数。

生存分析中，时间必须统一单位，事件 1=死亡，0=删失。随包真实数据用月。
第 09 个案例使用 data.frame(time,event) 与 Ridge Cox：

```r
library(survival)
fit_cox <- function(x,y) {
  z <- x; time <- y$time; event <- y$event
  coef(coxph(Surv(time,event)~ridge(z,theta=10,scale=FALSE)))
}
predict_cox <- function(model,x) as.numeric(x%*%model)
cindex <- function(y,prediction,domain) {
  time <- y$time; event <- y$event; risk <- prediction
  as.numeric(concordance(Surv(time,event)~risk,reverse=TRUE)$concordance)
}
```

`reverse=TRUE` 对应风险越高、事件越早的 Cox 风险方向。
如果直接传 Surv 对象，回调也必须改成接受 Surv，而不是继续读取 `$time`。
缺乏事件或可比较对的验证折无法估计 C-index，应提前设计合理折分，而不是事后忽略。
教学示例固定惩罚强度以便运行；正式研究中需要说明惩罚参数如何确定，以及是否使用内部尺度标准化。

## 11. pooled 指标和分域指标

单个队列的 pooled C-index 把这个队列的全部合格样本放在一起，包括可能混合的组织学亚型。
这与把所有 GEO 队列合在一起再算一个 C-index 不同，也与各亚型 C-index 的平均值不同。

pooled 指标既包括同亚型内部的可比较对，也包括不同亚型之间的可比较对。
如果亚型本身与预后相关，模型可能从亚型差异获得 pooled 区分能力。
若希望说明同一种癌症内也预测得好，应另外报告各亚型指标。
各域指标等权平均是 macro，不能与 pooled 混称。

当前论文柱状图是每个保留外部队列的 pooled 结果；历史四张 λ 曲线保留了当时三折 macro 选择过程。
不能把这些曲线重新命名为 pooled，也不能默认为当时用 pooled 重新选择过 λ。

## 12. λ 诊断图怎么做

```r
diagnostic <- lambda_diagnostics(project_fit,tr$x,tr$domain)
diagnostic$summary
head(diagnostic$coordinates)
source(system.file("examples","06_lambda_fixed_pca.R",package="CrossDomainAdjust"))
```

诊断函数先在 λ=0 的矩阵上做一次均值中心化 PCA，不额外做特征 Z-score；
此后所有 λ 的数据投到同一套 PCA 轴上。第 06 个案例把所有面板坐标范围也设成一致。
返回值包含每个样本的 PC1/PC2 坐标、域标签、λ、PCA 基底、基线解释率和汇总指标。

域 R² 用完整特征空间的“域间平方和 / 总平方和”，域间项按样本量加权。
它不是前两个 PC 的 R²，也不是每个基因先标准化后 R² 的平均值。
本函数的中心距离是欧氏距离；旧 `diagnose_cross_domain()` 则报告除以特征数平方根后的距离，数值定义不同。

λ=1 时域均值 R² 接近 0，只说明拟合域的均值差异被消除。平移后协方差完全不变，
投影后也可能保留其他分布差异。外部数据采用冻结转换，甚至不保证均值完全对齐。
因此附图应该写“逐步削弱域均值差异”，不能写“所有域间差异都已消除”。
最优预后 λ 必须由验证决定，不能由这张几何图决定。

## 13. 两个真实案例

肺癌数据：TCGA LUAD/LUSC 共 996 个开发样本，GSE3141 共 110 个外部样本。
胶质瘤数据：TCGA LGG/GBM 共 664 个开发样本，GSE43378 共 50 个外部样本。
每个案例都固定使用既有研究中的 60 基因集合。来源和处理说明保存在 extdata/README.md。

```r
lung <- readRDS(system.file("extdata","lung_worked_data.rds",package="CrossDomainAdjust"))
dim(lung$train$x)
table(lung$train$domain)
head(lung$train$survival)
source(system.file("examples","14_real_lung_survival.R",package="CrossDomainAdjust"))
source(system.file("examples","15_real_glioma_survival.R",package="CrossDomainAdjust"))
```

这两个脚本比较原始表达和四种方法，拟合 Ridge Cox，并输出外部 pooled C-index。
为了突出包的使用，示例固定 λ=.4、θ=10；它们不重跑基因发现，也不复现论文筛选后冻结模型的具体数值。
不要把示例输出复制到论文当成原始冻结结果。代码没有要求某个方法必须最好。
这些公开队列此前已被研究过程查看过，不能因为重新跑了例子就称为全新确认性验证。

换成自己的数据时，需要保留 x、domain、survival 的组织方式，核对基因、标签和终点。
如果外部平台缺少模型基因，应在新的开发分析中重新确定面板并重新拟合，而不是悄悄删基因继续预测。

## 14. 多个种子与结果保存

第 11 个案例用 11、22、33 三个固定种子，为中心平移和投影保留全部六个 λ 的得分。
同时标记最优 λ 是否在开区间 (0,1)，但不会删除端点结果。
正式分析可以预先增加重复次数，查看 λ 选择和性能是否稳定。

如果只展示两个 λ 同时非 0/1 的种子，需要在文章中说明这个筛选过程，并提供所有种子和完整网格。
反复找外部队列直到某种方法第一属于探索性筛选，不能按未经筛选的外部验证来解释。
保留全部候选队列及排除理由，有利于审稿时说明结果边界。

## 15. 保存模型和单样本部署

完整模型应保存预处理对象、下游模型、基因顺序、域定义、λ、惩罚参数、表达尺度和软件版本。
生存模型还需要说明事件方向和时间单位。

```r
prep <- project_fit
z <- transform_cross_domain(prep,tr$x)
model <- fit_model(z,tr$regression)
bundle <- list(preprocessing=prep,model=model,genes=colnames(tr$x),
               version=as.character(packageVersion("CrossDomainAdjust")))
path <- tempfile(fileext=".rds")
saveRDS(bundle,path)
saved <- readRDS(path)
new_z <- transform_cross_domain(saved$preprocessing,te$x[1,])
new_prediction <- predict_model(saved$model,new_z)
unlink(path)
```

变换时可以传 `lambda=` 临时覆盖拟合对象的强度，原对象不会被修改，适合画诊断图。
但不能在部署时随意改 λ，却继续沿用另一个 λ 下拟合的预测系数。
改变转换规则后，下游模型也需要在相应训练表示上重新拟合、验证。

## 16. 高级接口

`fit_cross_domain()` 允许分别控制 representation、standardize、correction。
顺序始终是：编码 → 可选标准化 → 域校正。其默认标准化是 none，而 rank/pair 预设显式启用 zscore。
预设固定的参数不能通过 `...` 覆盖；需要自定义时使用通用接口。

```r
hybrid <- fit_cross_domain(tr$x,tr$domain,representation="raw",standardize="none",
  correction="hybrid",lambda=c(.2,.4))
z <- transform_cross_domain(hybrid,te$x,te$domain)
```

hybrid 先平移再投影，两个 λ 按平移、投影排列；标量则两步共用。
这两步有重叠，不假设组合一定更好。它是高级接口，不是论文四种主要方法之外新增的必选模型。

底层函数还有 fit_feature_encoder/encode_features、fit_feature_scaler/scale_features、
fit_domain_adjuster/adjust_features，以及保留兼容的投影函数名。
底层单样本投影封装返回命名向量；完整 pipeline 的 transform 返回数值矩阵。

## 17. 常见问题

| 问题 | 处理方式 |
|---|---|
| 提示缺少基因名 | 确认列为基因并赋予唯一名称；很多数据需要转置 |
| 提示缺少拟合基因 | 提供完整冻结集合，或在开发阶段重新设计模型 |
| 输入有 NA/Inf | 在上游按训练规则解决；包不会自动填补 |
| 配对数量超上限 | 缩小训练基因面板或显式候选对，不能只设 max_pairs |
| 没有配对被保留 | 查看训练频率和候选定义，不在测试集上放宽筛选 |
| 未知域无法平移 | 需要已拟合域；投影可计算，但不保证去除新域差异 |
| 排序后超出 0–1 | 预设又做了冻结 Z-score，这是正常的 |
| λ 越大预测越差 | 可能去掉了预后相关信号，按验证结果选择 |
| CV 分数不是有限数 | 检查结局对应、事件数和可比较样本对，不应忽略失败折 |
| 外部 λ=1 仍有差异 | 冻结方法没有重新估计外部均值，不保证全部校正 |

2.1.0 保留 2.0.1 冻结转换行为；新增预设、诊断和案例不会改写旧论文结果。
2.0.0 没有 scaler 的对象按不额外标准化处理。0.1.x 旧序列化投影对象建议重新拟合。
历史源码包和冻结模型应保留，避免为了升级覆盖可复现依据。

## 18. 一次运行全部案例

```r
source(system.file("examples","run_all.R",package="CrossDomainAdjust"))
```

01 四方法；02 排序与缩放；03 配对筛选；04 已知域平移；05 单样本投影；
06 λ 图；07 回归 CV；08 分类；09 生存 CV；10 多域稳健中心；11 全种子网格；
12 保存与基因校验；13 高级组合；14 真实肺癌；15 真实胶质瘤。

前 13 个案例中的数值来自明确标注的模拟教学数据，不用于临床结论。
真实数据案例用于演示可执行流程。正式研究请保存来源、样本纳排、基因映射、预处理、完整种子网格、
结局定义、各队列选择记录、软件版本及最终冻结对象，确保研究结果可以逐项核对。
