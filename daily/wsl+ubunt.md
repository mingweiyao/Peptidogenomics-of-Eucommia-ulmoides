# Windows 上使用 WSL 安装 Ubuntu（可指定安装路径）

本文记录在 Windows 上使用 **WSL 2** 安装 **Ubuntu 24.04** 的完整流程，  
**重点：支持自定义安装路径（非 C 盘）**，适合开发 / 网络 / 实验环境。

---

## 一、环境说明

- Windows 11
- WSL 2
- Ubuntu 24.04
- 安装方式：`wsl --install`（可控路径）

---

## 二、安装 Ubuntu 并迁移路径

### 1. 以管理员身份打开 PowerShell

执行：

```powershell
# 查看可用发行版
wsl --list --online

# 安装 Ubuntu 24.04（会默认装到 C 盘）
wsl --install -d Ubuntu-24.04

# 备份为 tar
wsl --export Ubuntu-24.04 D:\ubuntu\Ubuntu2404_backup.tar

# 注销原实例
wsl --unregister Ubuntu-24.04

# 重新导入到指定路径（WSL 2）
wsl --import Ubuntu-24.04 D:\ubuntu\Ubuntu2404 D:\ubuntu\Ubuntu2404_backup.tar --version 2
```

### 2. 配置普通用户
```
wsl -d Ubuntu-24.04 -u root
adduser ymw
usermod -aG sudo ymw
ubuntu2404 config --default-user ymw
```

---

## 三、使用技巧

### 1. 切换到windows路径
a. 根目录在/mnt下，因此切换路径时要先处理/mnt，其中 磁盘符小写
b. 示例：cd /mnt/d/ubuntu/tools

### 2. conda 的安装

```powershell
# 删除安装目录及配置文件
rm -rf ~/anaconda3 ~/.condarc ~/.conda ~/.continuum

# 更新系统并安装必要依赖
sudo apt update && sudo apt upgrade -y
sudo apt install -y bzip2 wget

# 确认架构
uname -p

# 下载安装包并安装
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-2025.12-2-Linux-x86_64.sh
chmod +x Anaconda3-2025.12-2-Linux-x86_64.sh
bash Anaconda3-2025.12-2-Linux-x86_64.sh

# 激活配置
source ~/.bashrc
```

### 3. 配置镜像
```powershell
# 配置国内镜像源，加速包下载
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --set show_channel_urls yes
```

### 4. conda 安装和卸载包
```powershell
conda install -c bioconda sra-tools
conda uninstall sra-tools

# 更新或者重装包
conda install -c conda-forge -c bioconda --force-reinstall "sra-tools>=3.2.0"
```

### 5. conda 虚拟环境
```powershell
# 创建虚拟环境
conda create -n rnaseq

# 激活和关闭虚拟环境
conda activate rnaseq
conda deactivate

# 查看创建的虚拟环境
conda env list

# 删除指定环境
conda remove -n rnaseq --all
```

### 6. 挂载硬盘
```powershell
sudo mkdir -p /mnt/f
sudo mount -t drvfs F: /mnt/f
```
---

## 四、软件安装

### 1. aspera
```powershell
# 安装ascp
wget https://download.asperasoft.com/download/sw/connect/3.9.1/ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
tar zxvf ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.sh
echo 'export PATH=$PATH:/home/ymw/.aspera/connect/bin' >> ~/.bashrc
source ~/.bashrc

# 使用
tail -n +2 data_list.tsv \ # data_list.tsv 从ENA下载
| while IFS=$'\t' read -r run fastq_aspera sra_aspera; do
    echo "Downloading $run ..."
    ascp -P 33001 -k 1 -l 100m -T \
      -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
      "era-fasp@${fastq_aspera}" ./
  done

```

### 2. fastx-toolkit
```powershell
# 下载及解压安装
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
tar -xjvf fastx_toolkit-0.0.14.tar.bz2
cd fastx_toolkit-0.0.14
./configure
make

# 报错
make clean 2>/dev/null || true
FLAGS="-O2 -g -Wno-implicit-fallthrough" CFLAGS="-O2 -g -Wno-implicit-fallthrough" ./configure
sudo make install
```