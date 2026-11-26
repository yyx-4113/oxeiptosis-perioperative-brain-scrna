# GitHub 归档操作手册（中文）

## 1. 仓库信息
- 仓库名：`yyx-4113/oxeiptosis-perioperative-brain-scrna`（kebab-case）
- 默认分支：`main`
- 发布标签：`v1.0.0`

## 2. 首次创建
1. 在 GitHub 网页端以 README 初始化新建空仓库。
2. 本地进入 `analysis/repro_repo`，执行：
   - `git init -b main`
   - `git remote add origin https://github.com/yyx-4113/oxeiptosis-perioperative-brain-scrna.git`
   - 将 `../scripts/`、`../00_environment/`、`../01_data/meta/` 之外的产物目录按需复制进仓库根目录（保留 README/CITATION/LICENSE 的相对路径结构）。
3. 确认 `.workbuddy/` 与原始 GEO TAR 不入库。

## 3. 发布
- `git add -A && git commit -m "Release v1.0.0: two-cohort oxeiptosis reanalysis"`
- `git tag -a v1.0.0 -m "Reproduce manuscript v1.0"`
- `git push origin main --tags`
- 在 GitHub Releases 页面基于该 tag 创建 release。

## 4. 校验
- 下载 release 后运行 `python scripts/50_manifest.py` 与 `python scripts/90_audit.py`，确认 `audit_report.json` 中 `n_fail = 0`。
- 数据可用性声明必须写实名仓库 URL + 版本 tag，禁止 "available on request"。
