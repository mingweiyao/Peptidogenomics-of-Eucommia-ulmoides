import os, sys, platform, subprocess, urllib.request, tarfile, shutil

ROOT = os.path.abspath(os.path.dirname(__file__))
BIN_DIR = os.path.join(ROOT, ".mamba")
ENV_DIR = os.path.join(ROOT, ".env")
MAMBA = os.path.join(BIN_DIR, "micromamba.exe" if os.name == "nt" else "micromamba")
ENV_YML = os.path.join(ROOT, "environment.yml")
MAIN_PY = os.path.join(ROOT, "workflow_family_finding.py")

def run(cmd, **kw):
    print(">", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)

def download_micromamba():
    os.makedirs(BIN_DIR, exist_ok=True)

    if os.name == "nt":
        url = "https://micro.mamba.pm/api/micromamba/win-64/latest"
        archive = os.path.join(BIN_DIR, "micromamba.tar.bz2")
        urllib.request.urlretrieve(url, archive)
        with tarfile.open(archive, "r:bz2") as tf:
            tf.extractall(BIN_DIR)
        found = None
        for r, _, files in os.walk(BIN_DIR):
            for f in files:
                if f.lower() == "micromamba.exe":
                    found = os.path.join(r, f)
                    break
            if found:
                break
        if not found:
            raise RuntimeError("micromamba.exe not found after extraction")
        shutil.copyfile(found, MAMBA)

    else:
        # linux-64 / osx-64 / osx-arm64
        sysname = platform.system().lower()
        machine = platform.machine().lower()
        if sysname == "darwin" and machine in ("arm64", "aarch64"):
            plat = "osx-arm64"
        elif sysname == "darwin":
            plat = "osx-64"
        elif sysname == "linux":
            plat = "linux-64"
        else:
            raise RuntimeError(f"Unsupported platform: {sysname} {machine}")

        url = f"https://micro.mamba.pm/api/micromamba/{plat}/latest"
        archive = os.path.join(BIN_DIR, "micromamba.tar.bz2")
        urllib.request.urlretrieve(url, archive)
        with tarfile.open(archive, "r:bz2") as tf:
            tf.extractall(BIN_DIR)
        found = None
        for r, _, files in os.walk(BIN_DIR):
            for f in files:
                if f == "micromamba":
                    found = os.path.join(r, f)
                    break
            if found:
                break
        if not found:
            raise RuntimeError("micromamba not found after extraction")
        shutil.copyfile(found, MAMBA)
        os.chmod(MAMBA, 0o755)

def ensure_env():
    env = os.environ.copy()
    env["MAMBA_ROOT_PREFIX"] = ENV_DIR
    # 先创建（如果已存在会失败）
    try:
        run([MAMBA, "create", "-y", "-f", ENV_YML, "-n", "genefam"], env=env)
    except subprocess.CalledProcessError:
        # 已存在则更新
        run([MAMBA, "env", "update", "-y", "-f", ENV_YML, "-n", "genefam"], env=env)


def run_pipeline():
    env = os.environ.copy()
    env["MAMBA_ROOT_PREFIX"] = ENV_DIR

    # 用 micromamba run 在环境里执行 python run.py
    run([MAMBA, "run", "-n", "genefam", "python", MAIN_PY], env=env)

if __name__ == "__main__":
    if not os.path.exists(ENV_YML):
        raise SystemExit("environment.yml not found.")
    if not os.path.exists(MAIN_PY):
        raise SystemExit("workflow_Family_finding.py not found. (Set MAIN_PY ...)")
    ensure_env()
    run_pipeline()
