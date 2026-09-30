# Code availability

## Manuscript-ready statement

All code used to reproduce the statistical analyses and figures, together with
the frozen analysis definitions, analysis-ready derived tables, figure source
data and verification utilities, is publicly available at
https://github.com/FENG1567/acute-pancreatitis-ards-microbiome (release
`v1.0.1`) under the MIT License. The default workflow reproduces the primary
and secondary statistical results without downloading raw sequencing reads.
Raw sequencing and transcriptomic data are not redistributed in this
repository; they remain available from the original public repositories under
the accessions listed in `docs/ACCESSIONS.md`.

## Reproduction entry point

```bash
python -m pip install -r requirements.txt
python run_all.py --no-figures
python -m unittest discover -s tests -v
```

Run `python run_all.py --require-figures` when R and the required graphics
packages are installed and figure regeneration is required.

## Scope and licensing

The MIT License applies to the repository code. Public-source and derived
scientific data retain the terms of their originating repositories and
publications, as described in `DATA_USE_NOTICE.md`.

## 中文核对

用于复现统计分析和图片的代码、冻结分析定义、分析就绪的派生表格、图片源数据及验证工具，均已在
https://github.com/FENG1567/acute-pancreatitis-ards-microbiome 公开发布，版本为
`v1.0.1`，代码采用 MIT License。默认流程无需重新下载原始测序数据，即可复现主要及次要统计结果。
原始测序和转录组数据未在本仓库中重复分发，可根据 `docs/ACCESSIONS.md` 中列出的登录号从原始公共数据库获取。
