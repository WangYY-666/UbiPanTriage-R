# 网页端预计算文库（Web lookup library）

网页版页面 2-4（泛素单基因 / 自定义基因集 / 自定义单基因）不再加载整包癌种缓存，
而是读取"按基因分文件"的预计算文库，单次查询只读 1 个小文件，2 核 4G 服务器可秒级返回。

## 产物位置

`inst/extdata/lookup/`

- `gene/<GENE>.rds`  每个基因一个文件：`per_cancer[[癌种]]` = 评分行（含 Basic/Immune/
  Metabolic 及各分项、`Imm_<方法>_<特征>`、`Meta_<通路>` 相关列）+ 表达切片
  （TPM / count / FPKM，int16 编码 `round(log2(x+1)*100)`，-32768 = NA）
- `genes_index.rds`   全部基因清单（Gene + 覆盖癌种，压缩后 <1 MB）
- `sample_info.rds`   各癌种样本与肿瘤/正常标签
- `clinical.rds`      各癌种生存数据（sample / OS / OS.time）
- `.partial/<癌种>.rds` 中间结果，可复用（存在则跳过该癌种，`--rebuild` 强制重算）

## 构建命令

在项目根目录（含 `pancaner_data/`）：

    # 全量 33 癌种，16 核（本机内存需 >= 16G）
    UBI_CORES=16 Rscript UbiPanTriage/data-raw/build_web_lookup.R
    # 只构建部分癌种
    Rscript UbiPanTriage/data-raw/build_web_lookup.R --only=LUAD,BRCA
    # 小规模冒烟（前 800 个表达基因）
    Rscript UbiPanTriage/data-raw/build_web_lookup.R --only=ACC --max-genes=800

参数：`--min-tpm`（默认 1）、`--min-frac`（默认 0.2，即 >=20% 肿瘤样本表达）。
过滤后每癌种约 1.8-2.1 万个表达基因；全部 33 癌种构建约需 2-3 小时。

## 网页端运行数据

`app.R` 通过 `options(ubi.data_dir)` 或包内 `inst/extdata` 定位数据。运行时只需要：

- `pancancer_ubi_scores.rds`（页面 1-2 的泛素评分表，53 MB）
- `lookup/`（页面 2-4 文库，全 33 癌种约 1-2 GB）

旧的 `*_cache.rds`（共 2.2 GB）网页端不再使用，仅 R 包完整功能（实时重算、权重微调、
离线批量）需要，部署时可不随网页端上传。

## 数值约定（与 R 包一致）

- 表达切片为 int16 编码，包内 `gene_expr_raw()` 按 `2^(v/100)-1` 解码。
- 免疫列名 `Imm_<方法>_<特征>`（方法小写、无 `_TIMER` 等法后缀）；代谢列 `Meta_<通路>`。
- 基础分 = 表达/差异/生存/ROC 四个分项的加权均值（网页端默认等权）；
  无正常组织或差异数据的癌种在可用分项上重归一化并在 `Missing` 列标记。
