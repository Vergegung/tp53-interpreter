# TP53 Variant Atlas v0.4

TP53 结构、功能与显性负性的分层研究评估。浏览器直接计算，GitHub Pages 托管，无需用户启动 Python、安装软件或运行本地服务。

新版统一显示：原始值与来源内百分位 → 按实验家族汇总的结构 / 功能 / DNE 分项 → 可调整权重的综合研究指数。缺失与多记录歧义传播为区间，非 missense 进入序列机制路线。详细公式、来源分母、适用性与局限见 [METHODS.md](docs/METHODS.md)。

综合指数是探索性汇总，尚未外部校准；它不是致病概率、临床风险或药物响应概率。PDO、RNAseq 与临床验证后置。

## 在 GitHub 执行与发布

1. 将本仓库文件上传到 GitHub，保留 `docs/`、`tests/`、`.github/workflows/pages.yml`。
2. 仓库 Settings → Pages → Source 选择 **GitHub Actions**。
3. Push 到 main/master 或在 Actions 手工运行 **Validate and publish TP53 Atlas**。
4. 校验通过后发布；在 Actions 的 deployment 输出或 Settings → Pages 打开网站。

GitHub Pages 用户只访问网站，无需本地执行。测试使用 Node 原生 API，无 npm 安装。GitHub 官方说明：[Pages](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages)、[Actions 工作流](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages)。

## 数据和复现

`docs/data/public.json.gz` 仅包含公开分子数据：7,467 个错义替换的原始实验值 / 分数、23 个替换的结构实测、公共序列与外显子、来源分类与引用。患者临床、PDO 及 RNAseq 表未读取或导出。访问时解压并核对 `manifest.json` 的 SHA256，查询与手工记录仅在浏览器内处理。

来源继承已核验的 v0.3 快照，网页计算与原 Python 核心通过全部 SNV / 单碱基 indel 和代表轴的比对。代码升级与来源快照升级分别记录；数据许可与引用仍归上游，本项目不为上游数据库重新赋予许可证。

```bash
node --test tests/*.test.mjs
```

`scripts/export_public.py` 是数据导出器，依赖相邻的已锁定 interpreter 开发工程，用于维护者刷新数据，访问者无需运行。测试夹具来自 v0.3 独立审计，仅含公共变异后果。

## 核心限制

- S 当前汇总 apo 与锌的来源内排名，DNA 原始结合值单列；不同实验尚未跨方法校准。
- F 三个家族等权，D 与 Giacomelli LOF 相关；预设权重必须进行敏感性评估。
- 非 missense 不套用错义实验分值；复杂 HGVS 与实际剪接产物仍可能未决。
- UMD/p53.fr 独立快照未导入；没有患者结局模型。
