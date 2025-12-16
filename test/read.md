GeneFamilyFinder/                # 项目根目录（随便叫什么）
├── bootstrap.py                 # 自举启动器（自动装环境 + 跑主程序）
├── environment.yml              # 依赖清单（blast / hmmer / muscle / biopython）
├── workflow_Family_finding.py   # 你的主 pipeline 脚本（真正干活的）
│
├── input/                       # 输入数据（建议）
│   ├── query.fa                 # 初始家族序列
│   ├── genome.fa                # 基因组序列
│   └── proteins.fa              # 蛋白序列
│
├── output/                      # 所有结果输出
│   ├── total_hits.fasta
│   ├── blastp.tsv
│   ├── tblastn.tsv
│   └── hmm.tbl
│
├── .mamba/                      # bootstrap 自动生成（不要手动改）
│   ├── micromamba.exe           # Windows
│   └── micromamba               # Linux / macOS
│
├── .env/                        # micromamba root prefix（自动生成）
│   ├── envs/
│   │   └── genefam/             # 真正的 conda 环境
│   ├── pkgs/                    # 包缓存
│   └── condarc / metadata
│
├── run.bat                      # （可选）Windows 双击入口
├── run.sh                       # （可选）Linux/macOS 启动脚本
│
└── README.md                    # （强烈建议）一句话用法说明
