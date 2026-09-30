# 急性胰腺炎跨部位微生物与早期 ARDS：代码复现仓库

本仓库包含论文的冻结分析方案、分析就绪数据、统计代码、图源数据、最终 PDF 图和自动核验工具。

公开仓库地址为：<https://github.com/FENG1567/acute-pancreatitis-ards-microbiome>。
论文候选版本标记为 `v1.0.0`，代码采用 MIT License。可直接用于论文的中英文代码可用性说明见
`CODE_AVAILABILITY.md`。

最简复现命令：

```powershell
python -m pip install -r requirements.txt
python run_all.py
```

默认流程可从仓库内的小型计数矩阵重算 PRJNA893348 的 65 名 AP 患者主分析、敏感性分析、Firth 模型，以及 PRJNA428535 的 62 对血液–中性粒细胞分析。若本机有 R、`ragg`、`svglite` 与 `systemfonts`，还会重绘六张主图。

必须注意：主检验 P 值是 `0.009299070092990702`。`sensitivity_with_ci.tsv` 中的 `0.008799120087991202` 来自另一条固定随机种子敏感性随机流，不能替代主 P 值。

原始 FASTQ、DADA2 RDS、RDP 数据库和完整宿主表达矩阵体积较大，未放入 GitHub。复现者应按 `docs/ACCESSIONS.md` 从公共数据库重新下载。详细边界见 `docs/REPRODUCIBILITY_SCOPE.md`，本次整理与补齐情况见 `docs/GITHUB_PACKAGE_REPORT_zh-CN.md`。
