# 交接说明 · Week 5 生物信息学作业

> 给接手的助手:这份文档是自包含的。你看不到之前的对话,所有必要信息都在下面。
> 用户是 Windows 机器。请**先读完第 5 节的已知陷阱**再动手,那几个坑我已经踩过了。

---

## 1. 基本情况

| | |
|---|---|
| 学生 | Zhou Hanyi(周涵艺),学号 **SUAT24000202** |
| 课程 | Bioinformatics: From Multi-Omics Data to Discovery(SUAT 2026 秋) |
| 课程仓库 | `xielab2017/Bioinformatics_SUAT_2026_FALL`(公开) |
| 作业仓库 | `santianwan/Bioinformatics_homework_zhouhanyi`(学生自己的,分支 `main`) |
| 本机课程材料 | `F:\subject\3_up\bioinformation\Bioinformatics_SUAT_2026_FALL-main\` |
| 本机 EMP 工具 | `F:\EasyMultiProfiler-Web-9.0.4\EasyMultiProfiler-Web-9.0.4\` |
| 报告语言 | **英文**(与课程材料一致),对话用中文 |

Week 5 有两份作业。**HW1**(DESeq2 差异表达)与 **HW2**(在 EasyMultiProfiler-Web 网页端分析并用其 Sync 提交)。

---

## 2. 已经完成的部分

全部已推送到作业仓库 `Week5/` 目录(31 个文件)。结构:

```
Week5/
├── README.md                          总说明 + 提交清单勾选表
├── week5_interpretation.md            HW1 的 142 词规定解读 + 支撑细节
├── week5_AI_verification_log.md       HW1 的 AI 使用与验证日志
├── week5_homework2_guide.md           HW2 的操作指南 + 参考结果 + 解读草稿
├── code/
│   ├── week5_deseq2_analysis.R        ★ HW1 正式提交脚本(见下方警告)
│   ├── week5_truth_benchmark.R        HW1 事后校验(R 版,需先跑上面那个)
│   ├── week5_pydeseq2_crossrun.py     HW1 的 PyDESeq2 对照运行
│   ├── week5_figures_and_benchmark.py HW1 出图 + 事后校验
│   └── hw2_emp_rnaseq_reference.py    HW2 的独立参考分析
├── data/                              课程原始输入的副本(未修改)
├── figures/                           4 张图
└── outputs/                           结果表与摘要 JSON
```

Week 4 也已完成并提交在 `Week4/`,无需再动。

### ⚠ 最重要的一条:R 脚本从未被执行过

`week5_deseq2_analysis.R` 是正式提交物,内容完整正确,但**在我的沙箱里跑不起来**(原因见第 5.2 节)。为了让报告里的数字是算出来的而非编的,我用 **PyDESeq2 0.5.4**(DESeq2 方法的 Python 参考实现,含 apeglm 收缩)拟合了同一个设计。

**报告里现有的所有数字都来自 PyDESeq2。R 跑出来的结果为准。** 若两者有出入,请用 R 的结果更新 `week5_interpretation.md`,并在 `week5_AI_verification_log.md` §4 记下差异。DESeq2 与 PyDESeq2 高度一致但不逐位相同(Cook's 距离离群重拟合处尤其可能有差别)。

---

## 3. 剩下要做的两件事

### 任务 A:在本地跑 HW1 的 R 脚本,补齐三个缺失的提交物

作业要求 8 个文件,目前缺 3 个,它们只能由 R 脚本生成:
`week5_deseq2_results.csv`、`week5_deseq2_object.rds`、`session_info.txt`。

```powershell
cd D:\zhy
git clone https://github.com/santianwan/Bioinformatics_homework_zhouhanyi.git
cd Bioinformatics_homework_zhouhanyi\Week5\code

# 装包,只需一次
Rscript -e "if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org'); BiocManager::install(c('DESeq2','apeglm'),ask=FALSE,update=FALSE)"

Rscript week5_deseq2_analysis.R
```

脚本用 `../data/` 这样的相对路径,**工作目录必须是 `Week5\code`**。若用 RStudio,先
Session → Set Working Directory → To Source File Location,再 Source。

跑完会写出 `../outputs/` 的三个文件,并用 R 版本覆盖 `../figures/week5_pca.png` 和
`week5_de_plot.png`。然后:

```powershell
cd ..\..
git add Week5
git commit -m "Week 5 HW1: R outputs from local DESeq2 run"
git push
```

**跑出来应该对上的数字**(PyDESeq2 得到的,R 应当接近):

| 项目 | 值 |
|---|---|
| 过滤后参与检验的基因 | 989 / 1000(≥10 counts 于 ≥3 样本) |
| 显著(padj<0.05 且 \|log2FC\|≥1) | **60**(36 上调,24 下调) |
| 仅过 padj<0.05 | 83 |
| padj = NA | 0 |
| 最显著基因 | Gene0035,log2FC +1.74,padj 1.3×10⁻¹⁰ |
| PC1 / PC2 方差 | 24% / 9% |
| PC1 ~ condition | p = 4×10⁻¹² |
| PC1 ~ batch | p = 0.99(批次不在 PC1 上) |
| 系数名(R) | `condition_treated_vs_control` |

若 R 给出的显著基因数与 60 差得很远(比如相差一倍以上),说明某处设置不对,先查参考水平是不是 control、design 是不是 `~ batch + condition`。

### 任务 B:HW2,在 EasyMultiProfiler-Web 网页端完成并 Sync

**这一步必须由用户手动点击,任何助手都替不了**(原因见 5.3)。

要上传的两个文件(`tests/` 里只有这一套是 RNA-seq):

```
F:\EasyMultiProfiler-Web-9.0.4\EasyMultiProfiler-Web-9.0.4\tests\RNAseq_output.csv
F:\EasyMultiProfiler-Web-9.0.4\EasyMultiProfiler-Web-9.0.4\tests\RNAseq_mapping.csv
```

其余 `16S_*.csv`、`ChIP/*.bed`、`Clinical*.csv`、`meta-*.csv` 都不是 RNA-seq,不要传。

操作:双击 `Start-EMP-Web.bat` → 浏览器 `http://127.0.0.1:8080` → Analyze → Import 两个
CSV → Preprocess(记下过滤规则)→ 差异分析(参考水平 **DMSO**)→ PCA + 火山图 →
Export 页:学号 `SUAT24000202` + 姓名 + 口令(≥8 位)注册登录,需要时贴 GitHub PAT
(对作业仓库 Contents: Read and write)绑定 → 轨道 **transcriptomics**、作业 **Week 5** →
点同步 → 打开返回的 commit 链接核对。结束后双击 `Stop-EMP-Web-Windows.bat`
(关浏览器不停服务)。

结果会落在作业仓库的 `EMP2026/Week_05/transcriptomics/weekly/runs/<时间戳>/`,增量写入。

**完整的操作步骤、参考数字、解读草稿和同步前检查清单都在
`Week5/week5_homework2_guide.md`,直接读那份。**

HW2 数据集的关键点(容易做错的地方):`Group` 那一列看着是 6 个平行组,实际编码的是
**两个交叉因子**——3 化合物(DMSO 对照 / T4400 / T3976)× 2 超声(±LIPUS),每格 n=4。
参考结果:只有 T4400 有反应(145 个基因,82 上 63 下);T3976 和 LIPUS 都是 0;两个
交互项也都是 0(最小 padj 0.96 与 0.999),所以**没有证据支持超声增效**。

---

## 4. 两条必须守住的底线

**4.1 Instructor Key 的用法。** 课程包里的
`Week5_Homework_Gene_Annotation_Instructor_Key.csv` 含 `truth_log2FC_for_instructor` 列,
即模拟数据的真实效应量——这是答案,疑似老师误传进学生目录。

用户已明确同意的做法是:**分析全程不碰它,跑完之后只用于事后校验,并在报告里写明**。
现有代码严格遵守了这一点(`week5_deseq2_analysis.R` 和
`week5_pydeseq2_crossrun.py` 从不打开它,只有 `week5_figures_and_benchmark.py` 在分析
定稿后读取)。**请不要用它去挑阈值、筛基因或反推结果。**

已得到的校验结果:灵敏度 0.659,精确度 0.967,经验 FDR 0.033(名义 0.05),
58/58 方向全对,回归斜率 0.843。

**4.2 不要编造没跑过的结果。** 这份作业此前出过一次这类错误:在只检验了主效应的情况下
写了"交互不显著"的结论。后来补跑了交互检验,结论恰好成立,但**当时那句话是没有依据的**。
如果某个数字你没有实际算出来,就说没算,不要从相近的结果外推。

---

## 5. 已知陷阱(都已验证,别重复踩)

**5.1 F: 盘无法授予文件访问权。** 对 `F:\subject\...` 和
`F:\EasyMultiProfiler-Web-9.0.4\...` 的父目录与子目录各试过,全部返回
"cannot place its delete-withholding grant"。这是整块盘的问题,换文件夹没用。
**C: 和 D: 正常。** 需要读 F: 上的东西,请用户复制到 D: 再给路径。

好消息:课程材料和 EMP 的 `tests/` 数据**都在公开 GitHub 仓库里**,可以直接取,
不需要本机文件。我核对过课程 Week 5 的 7 个文件与仓库 blob SHA-1 完全一致。

**5.2 DESeq2 在受限沙箱里装不上(如果你也在沙箱里)。**
- bioconda **没有 win-64 构建**的 `bioconductor-deseq2`,conda 会静默跳过;
- 从 Bioconductor 下 Windows 二进制包能装进可写目录,但加载时报
  `Refusing to dyn.load shared library from writable path`——编译的 `.dll` 不允许从可写
  目录加载,而只读的 conda 库只能用 conda 改。

替代方案就是 PyDESeq2(`pip install pydeseq2`)。它在沙箱里还需要一个 joblib 补丁:
PyDESeq2 用 `inner_max_num_threads=1` 打开并行上下文,只有进程型后端接受该参数,而
loky 在沙箱里建不了命名管道。解法(已写进现有脚本):

```python
from joblib import register_parallel_backend
from joblib._parallel_backends import SequentialBackend
class _SB(SequentialBackend):
    supports_inner_max_num_threads = True
register_parallel_backend("sandbox_sequential", _SB)
INFERENCE = DefaultInference(n_cpus=1, backend="sandbox_sequential")
```

**用户本机不受这些限制**,本机装 DESeq2 是正常的,任务 A 直接跑就行。

**5.3 `http://127.0.0.1:8080` 助手访问不了。** 私有/保留 IP 是硬性安全拒绝,连申请授权
的入口都没有;而且沙箱自己的回环地址不是用户的机器。EMP 的上传、网页分析、Sync 按钮
**只能用户手动操作**。不要试图绕过,也不要假装跑过网页端的分析。

**5.4 其他小坑。** PowerShell 沙箱里没有 `gh` 也没有 `python` 命令(用 GitHub REST API
和 `sys.executable`);`api.anaconda.org` 在拒绝名单上查不了包;Bioconductor 的下载会重定向到
CDN `mghp.osn.xsede.org`,需要单独放行。

---

## 6. 一句话总结

代码、图、解读、验证日志都做好并已推送;**剩下两件事都必须在用户本机执行**:
跑一次 `Rscript week5_deseq2_analysis.R` 补齐三个 R 输出文件,以及在
EasyMultiProfiler-Web 网页端完成 HW2 并点 Sync。先读
`Week5/README.md` 和 `Week5/week5_homework2_guide.md`。
