# 浏览器验收 v0.4.0

日期：2026-10-06。浏览器：Codex 内置浏览器；本地 HTTP 仅用于开发验收，最终目标为 GitHub HTTPS 静态站点。

- 数据加载：7,467 条蛋白替换，gzip 解压后的 SHA256 与 manifest 一致后才开启解读。
- 默认：R175H 76.4–80.9、R273H 69.8–71.9、Y220C 72.8–86.1；覆盖分别 100%、100%、87.5%（界面舍入 88%）。
- 三项等权：R175H 76.0–79.0、R273H 63.6–65.0、Y220C 69.6–86.8；原始值和分层值保留。
- 来源表：单位、原始范围、方向统一后的百分位和 N=23/20/2314/7448/7448/7448/4069 均显示。
- 非 missense 示例：c.586C>T、c.563del、c.376-2A>T 均“尚未定量”；不补 0、不推定零蛋白。无义记录保留 PTC 外显子 6、位置性 NMD 以及等位基因/剪接门控。
- 补充测量：多个变异时拒绝加入；单个 R175H 的模拟 Tm=36°C 单列为未核验，指数仍为 76.4–80.9。模拟记录已移除。
- 参考残基错误 R176H：明确报错，原结果隐藏，避免旧结果被误认为新查询。
- JSON、SVG、Markdown 下载：内置浏览器下载事件接口未返回事件，但三个实际下载文件已在 Downloads 核对并复制到 browser_exports/；JSON 数值与引擎一致，SVG 可解析，Markdown 包含评分和限制。
- 桌面：内容宽 1265 px，无页面水平溢出；卡片、权重、来源与机制三个视图正常。截图 TP53_unified_v04_desktop.jpg。
- 窄屏：390×844 px 视口，内容宽 375 px（扣除滚动条），页面滚动宽也是 375 px；输入、按钮、比较卡片改为纵向。截图 TP53_unified_v04_mobile.jpg。
- 浏览器 console error/warn：空。
- 临时视口已恢复；本地预览服务仅验收使用。

2026-10-06 已补充线上验收：[https://vergegung.github.io/tp53-interpreter/](https://vergegung.github.io/tp53-interpreter/)。GitHub Actions validate/deploy均成功，云端13组测试通过；Safari实际加载、SHA256数据门控、默认三突变评分、三项等权切换、非错义机制和三种实际下载文件核对通过。在线浏览器 error/warn 为空。部署和云端测试日志另存于 qa/cloud_*_20261006.txt；线上截图与下载记录位于本地 qa/deployment_20261006/。

云端记录：[https://github.com/Vergegung/tp53-interpreter/actions/runs/37434429493](https://github.com/Vergegung/tp53-interpreter/actions/runs/37434429493)。
