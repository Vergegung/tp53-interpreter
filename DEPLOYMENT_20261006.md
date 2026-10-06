# TP53 解读器 v0.4 已上线

在线网站：[https://vergegung.github.io/tp53-interpreter/](https://vergegung.github.io/tp53-interpreter/)

源代码：[Vergegung/tp53-interpreter](https://github.com/Vergegung/tp53-interpreter)

[GitHub 云端测试与发布记录](https://github.com/Vergegung/tp53-interpreter/actions/runs/37434429493)：validate 与 deploy 均成功，13 组测试全部通过。23 个初始上传文件的 Git blob SHA 与已核验本地内容一致。

线上数据校验、默认 R175H/R273H/Y220C 评分、三项等权切换、非 missense 机制、JSON/SVG/Markdown 实际下载均已验证。Safari 可正常加载；后台浏览器验收 error/warn 为空。

访问网站即可使用，无需 Python 或本地服务器。GitHub Actions 负责测试和发布，评分在访问者浏览器内完成。患者 PDO、RNAseq 与临床数据未上传；指数仍为未经临床外部校准的研究评估。

原始发布提交：f2e6eae8382d7fcf2c80c809c75fac9d9bcbc559。分子来源快照：df49d0a729ef02bb7661df2bce7b9245fdfbf6b77c88017615a6d324b16b676f。

原来的桌面 v0.4 归档保留发布前状态；本记录更新真实上线状态，不改变其原始校验和。
