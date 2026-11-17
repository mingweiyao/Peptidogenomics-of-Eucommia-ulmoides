from pathlib import Path
import pymupdf, os, io
from PIL import Image

B = str(Path(__file__).resolve().parents[2] / "output" / "04_article")
JOBS = {
    "R3_transcript_architecture/Fig4": ["R3_transcript_architecture/figures/Fig4A-G.pdf",
                                        "R3_transcript_architecture/figures/Fig4H-I.pdf"],
    "R4_tissue_specificity/Fig5": ["R4_tissue_specificity/figures/Fig5A-E.pdf",
                                   "R4_tissue_specificity/figures/Fig5F-G.pdf",
                                   "R4_tissue_specificity/figures/Fig5H_examples.pdf"],

}
GAP = 12   

for out_name, parts in JOBS.items():
    docs = [pymupdf.open(f"{B}/{p}") for p in parts]
    W = max(d[0].rect.width for d in docs) + 2
    H = sum(d[0].rect.height for d in docs) + GAP * (len(docs) - 1) + 2
    out = pymupdf.open()
    page = out.new_page(width=W, height=H)
    y = 1.0
    for d in docs:
        r = d[0].rect
        x = (W - r.width) / 2
        page.show_pdf_page(pymupdf.Rect(x, y, x + r.width, y + r.height), d, 0)
        y += r.height + GAP
        d.close()
    out_path = f"{B}/{out_name.split('/')[0]}/figures/{out_name.split('/')[1]}.pdf"
    out.save(out_path)
    pix = page.get_pixmap(matrix=pymupdf.Matrix(600/72, 600/72))
    img = Image.open(io.BytesIO(pix.tobytes("png")))
    tdir = f"{B}/{out_name.split('/')[0]}/figures/tiff_600"
    os.makedirs(tdir, exist_ok=True)
    img.save(f"{tdir}/{out_name.split('/')[1]}.tiff", compression="tiff_deflate")
    out.close()
    print(f"merged -> {out_path} ({W:.0f}x{H:.0f} pt) + tiff")
print("done")
