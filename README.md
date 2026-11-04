# Binary prediction modeling and variable importance

将二分类模型开发、独立时间验证、变量重要性提取及 Top10 绘图串联起来的 R 脚本项目。支持 GBM、GLM、RF、DL、xGboost、Bayes 和 SVM，共 7 种算法；内置 12 套临床及病毒特征方案。

仓库只包含代码、配置和模拟数据生成器。真实病例数据、已经训练好的模型及分析结果需要在本地提供，未包含在此文件夹中。

## 1. 环境与依赖

- R ≥ 4.1；`ggplot2` ≥ 3.4；需支持 Cairo 图形设备。
- R 包：`h2o`、`e1071`、`pROC`、`PRROC`、`ggplot2`、`stringr`。
- H2O 模型训练还需要兼容所用 H2O 版本的 Java；纯 SVM 示例不启动 H2O 服务。
- xGboost 使用 H2O 后端，需要 `h2o.xgboost.available()` 返回 TRUE；不支持时会报错，不会自动删去该算法。
- 加载旧 H2O 二进制模型时，R 客户端和 H2O 服务端版本应与保存模型时一致。安装脚本只补齐缺失依赖，不替换已经安装的版本。

在仓库根目录运行：

```bash
Rscript scripts/install_dependencies.R
```

本次整合验证所用版本见 [docs/validation.md](docs/validation.md)。这里未生成 `renv.lock`，因此没有承诺不同依赖版本之间结果完全一致。

## 2. 快速运行模拟示例

```bash
Rscript scripts/make_example_data.R
Rscript run_pipeline.R config/example.R check
Rscript run_pipeline.R config/example.R all
```

示例完全随机模拟，包含 80 例开发数据、28 例时间验证数据，与真实队列无关。为了快速验证，示例仅运行 1 套特征、径向核 SVM、3 折交叉验证、200 次性能 CI bootstrap 和 10 次重要性置换。正式配置的设置见下一节。

输出位于 `outputs/example/model_outputs/`，其中 `variable_importance/` 包含重要性表及 Top10 PDF、PNG。示例数据和输出均已加入 `.gitignore`。重复运行时请使用新的输出目录；脚本不会自动清除已有结果。

## 3. 使用自己的数据

将输入表放入本地 `data/private/`，然后修改 `config/project.R`。也可以复制为 `config/local.R` 后修改；该本地配置默认不进入 Git。

### 模型数据表

`data_file` 指向以制表符分隔的表格，一行对应一名患者，至少包含：

| 列 | 含义 |
|---|---|
| `sample_code` | 唯一患者 ID，由 `id_col` 指定 |
| `Group` | 二分类结局，由 `outcome_col` 指定；默认阳性标签为 `CIP` |
| 预测变量列 | 必须与所选方案的 `predictors` 完全对应 |

阴性标签由开发集的另一类自动确定。两组结局标签必须一致，必须各包含两个类别，开发集每类样本数至少为 `nfolds`。不需要预先把结局转换为 factor。

### 队列划分表

`labels_file` 至少包含唯一 ID 列和队列标签列。例如：

```text
sample_code	label
SIM_0001	train
SIM_0002	test
```

`label_id_col = NULL` 表示使用划分表第一列作为 ID；默认标签列为 `label`，开发集为 `train`、时间验证集为 `test`。两张表的 ID 集合必须完全一致，不能有重复或未分配患者。患者不能同时进入两个队列。

### 配置项

| 配置 | 用途 |
|---|---|
| `data_file` / `labels_file` | 模型数据与队列划分表路径 |
| `output_dir` | 本次运行输出目录 |
| `models_root` | 默认 NULL，即 `<output_dir>/model_outputs`；也可指定已有模型输出目录 |
| `models` | 参与建模的算法列表 |
| `schemes` | 特征方案；每个方案含 `predictors` 和 `categorical_predictors` |
| `nfolds` / `seed` | 默认单次共享分层 5 折，随机种子 2026 |
| `thresholds` | 默认 NULL，在开发集 OOF 上选 Youden 阈值；指定数值则使用预设阈值 |
| `bootstrap_B` | 默认 2000，用于固定预测的性能 CI；不代表重拟合模型的 bootstrap 内部验证 |
| `parameters` | 各算法的预设参数；默认没有自动调参 |
| `importance_models` | 默认 5 种 H2O 原生重要性模型及 SVM；不含 Bayes |
| `importance_dataset` | SVM 置换数据，默认 `temporal`；可设 `development` |
| `n_perm` | SVM 每个原始变量的置换次数，默认 100 |

配置中的相对路径统一相对于仓库根目录，因此也可从其他工作目录调用脚本。

完整运行：

```bash
Rscript run_pipeline.R config/project.R check
Rscript run_pipeline.R config/project.R all
```

只运行部分方案时，例如在配置文件末尾添加：

```r
config$schemes <- config$schemes[c(
  "simplify_family_abun_model_01",
  "simplify_family_abun_model_02"
)]
```

## 4. 分阶段运行

```bash
# 开发、最终模型训练、独立时间验证
Rscript run_pipeline.R config/project.R model

# 从已有最终模型提取重要性，不重新训练
Rscript run_pipeline.R config/project.R importance

# 从已有重要性表生成独立 Top10 图，不重新计算分数
Rscript run_pipeline.R config/project.R plot
```

`all` 依次执行以上三个步骤。`check` 只验证输入数据和配置。

如果使用已有分析目录，在本地配置中设置 `config$models_root` 为该目录，并设置匹配的方案和算法，再运行 `importance` 或 `plot`。已有目录应保留完整 `development/manifest.rds`、`development/models/`、预处理信息以及预测表；仅有对外交付的模型文件夹不足以恢复原始变量映射。若重要性表已经存在，可直接运行 `plot`。

不应对已经完成的目录再次执行 `all` 或 `model`。提取和绘图也默认要求新目录；底层函数提供 `outdir` 参数，可在 R 中指定另一个输出目录。


### 加载保存的模型

以下命令在仓库根目录的 R 会话中运行。先根据配置定位已经完成的方案目录，并读取模型清单；清单会自动解析模型文件路径，兼容新运行的相对路径和已有运行的绝对路径。

```r
source("run_pipeline.R")
load_pipeline(normalizePath("."))

cfg <- resolve_config("config/project.R", normalizePath("."))
development_dir <- file.path(
  cfg$models_root, "simplify_family_abun_model_01", "development"
)
manifest <- read_model_manifest(development_dir)
```

**H2O 模型：使用 `h2o.loadModel()`。** GBM、GLM、RF、DL、xGboost、Bayes 均采用此方式；下面以 GBM 为例，将模型名称替换成对应算法即可。

```r
library(h2o)
h2o.init(max_mem_size = "4G", nthreads = 1L)

h2o_model <- h2o.loadModel(path = manifest$model_files[["GBM"]])
print(h2o_model)
```

**SVM 模型：使用 `readRDS()`。** 该模型由 `e1071` 训练，保存为 `SVM.rds`，加载时不需要启动 H2O 服务。

```r
library(e1071)

svm_model <- readRDS(manifest$model_files[["SVM"]])
print(svm_model)
```

如果只运行了模拟示例，将前面的配置文件改为 `config/example.R`，再执行 SVM 加载段；该示例没有训练 GBM。

加载后预测新数据时，配套预处理对象为 `manifest$preprocessor`：H2O 模型使用 `apply_native_preprocessor()`，SVM 使用 `apply_preprocessor()`，再通过 `predict_engine()` 获得阳性概率。预处理规则沿用训练时保存的对象。

## 5. 输出结构

```text
<output_dir>/model_outputs/
├── <scheme>/
│   ├── development/
│   │   ├── internal_oof/       # 每名开发患者在未用于拟合该折模型时的预测与评价
│   │   ├── apparent/           # 最终模型在全部开发数据上的表观评价
│   │   ├── models/             # 最终模型二进制文件
│   │   ├── manifest.rds        # 模型、配置、预处理与患者/折信息
│   │   └── ...                 # 汇总表、阈值、编码映射等
│   └── temporal/               # 锁定模型的独立时间验证
├── h2o_varimp/<scheme>/        # 原生重要性 TSV
├── svm_varimp/temporal/<scheme>/ # SVM 置换重要性 TSV；development 时目录同名替换
└── variable_importance/
    ├── README.txt
    └── <scheme>/
        ├── <model>_varimp.tsv  # 全部重要性条目
        ├── <model>_top10.pdf
        └── <model>_top10.png   # 300 dpi
```

原始评价输出每种模型包含指标表、预测表、阈值指标、ROC/PR/校准/决策曲线数据及独立 PDF、PNG。有关 CI 范围和解释限制，见 [docs/methods.md](docs/methods.md)。

重要性表包含 `variable`、`relative_importance`、`scaled_importance`、`percentage`。`percentage` 是 0—1 比例。原生 H2O 结果保留分类水平的独立行；SVM 按原始变量整体置换，保留负分。Top10 按标准化分数降序选取；不足 10 项则全部展示，负分为橙色。

## 6. 内置的 12 套特征方案

所有方案的精确变量列表由 `config/feature_schemes.R` 定义。

| 方案前缀 | 临床指标 |
|---|---|
| `simplify_*_model_01` | Smoking_history、Chronic_lung_disease_any、CD8_Ts、IL_8、autoAb_domain_positive |
| `simplify_*_model_02` | Smoking_history、Chronic_lung_disease_any、CD8_Ts、IL_2R、autoAb_domain_positive |
| `overall_*_model_01` | Smoking_history、Chronic_lung_disease_any、CD8_Ts、IL_8、CD4_Th、IL_2R、autoAb_domain_positive |

每类临床方案分别组合 `family_abun`、`family_det`、`species_abun`、`species_det`。`abun` 表示丰度，`det` 表示是否检出。

- 科水平：Herpesviridae。
- 种水平：HSV-1、EBV、HCMV、HHV-6A、HHV-6B、HHV-7，代码中的列名见特征配置。
- 吸烟史、慢性肺病、自身抗体分组及 `det` 指标按分类变量处理；其他指标按数值处理。

可以自行添加新方案，不要求运行恰好 12 套，也不要求队列恰好为 77/26 例。

## 7. 文件说明

| 文件 | 作用 |
|---|---|
| `run_pipeline.R` | 命令行统一入口，也提供 `load_pipeline()` |
| `config/project.R` | 正式运行配置，默认 7 种模型与 12 套方案 |
| `config/example.R` | 小规模模拟示例配置 |
| `config/feature_schemes.R` | 12 套方案的预测变量及分类变量列表 |
| `R/model_5fold_temporal.R` | 预处理、共享交叉验证、最终建模、时间验证与指标输出 |
| `R/importance_h2o.R` | H2O 原生重要性提取和变量名称还原 |
| `R/importance_svm.R` | 径向核 SVM 的 AUC 置换重要性 |
| `R/plot_importance.R` | 重要性表整理和独立 Top10 图 |
| `R/utils.R` | 数据检查、配置路径解析及模型恢复 |
| `R/pipeline.R` | 各阶段串联 |
| `scripts/make_example_data.R` | 生成纯模拟数据 |
| `scripts/install_dependencies.R` | 显式安装缺失的 R 包 |
| `tests/smoke_test.R` | 临时目录内的模拟全流程检查 |

## 8. 验证及 GitHub 上传

```bash
# 模拟数据下验证 SVM 全流程
Rscript tests/smoke_test.R

# 加上 H2O GBM 原生重要性分支，需要可用的 H2O 服务/Java
Rscript tests/smoke_test.R --h2o
```

将本文件夹中的内容作为仓库根目录上传即可，GitHub 会自动展示本 README。文件夹中未初始化 Git，也未连接或上传远程仓库。

`.gitignore` 排除了本地数据、运行结果、模型二进制、R 会话文件和本地配置；若通过网页上传，请只选择当前代码包中的文件，不额外拖入真实数据或生成结果。未替代码所有者选择开源许可证；公开发布前可由所有者补充 LICENSE。
