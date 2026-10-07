# import os
# os.environ['CUDA_VISIBLE_DEVICES'] = '' 
from ultralytics import YOLO
import time
import cv2
import base64
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr
from ultralytics.utils.plotting import colors 

def save_interactive_svg(result, svg_path, hover_line_width=3, font_size_ratio=0.02):
    """
    save interactive_svg
    """
    # 获取底图
    plotted = result.plot(labels=False, conf=False)
    h, w = plotted.shape[:2]

    # 动态计算文字大小：最短边的比例
    font_size = max(10, int(min(h, w) * font_size_ratio))

    ok, png_buf = cv2.imencode(".png", plotted)
    if not ok:
        raise RuntimeError("Failed to encode image to PNG.")
    b64 = base64.b64encode(png_buf.tobytes()).decode("ascii")
    data_uri = f"data:image/png;base64,{b64}"

    names = result.names if hasattr(result, "names") else {}
    boxes = []

    if result.boxes is not None and len(result.boxes) > 0:
        xyxy = result.boxes.xyxy.cpu().numpy()
        cls = result.boxes.cls.cpu().numpy().astype(int)
        conf = result.boxes.conf.cpu().numpy()

        for i in range(len(xyxy)):
            x1, y1, x2, y2 = map(float, xyxy[i])
            c = int(cls[i])
            boxes.append({
                "id": i,                                    # <-- 新增 box ID
                "x1": x1, "y1": y1, "x2": x2, "y2": y2,
                "class_name": names.get(c, str(c)),
                "confidence": float(conf[i]),
                "class_id": c
            })

    # 保存置信度到 txt，加一列 id
    with open(f"{svg_path}.confidence.txt", "w", encoding="utf-8") as f:
        f.write("id\tclass\tconfidence\n")                   # <-- 表头
        for box in boxes:
            f.write(f"{box['id']}\t{box['class_name']}\t{box['confidence']:.4f}\n")

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">')
    svg.append(f'  <image href={quoteattr(data_uri)} x="0" y="0" width="{w}" height="{h}" preserveAspectRatio="none" />')
    
    svg.append('  <style>')
    svg.append(f'    .det-box {{ fill: transparent; stroke: transparent; stroke-width: {hover_line_width}; cursor: pointer; }}')
    svg.append(f'    .det-group:hover .det-box {{ stroke: rgba(255, 255, 255, 0.8); }}')
    svg.append('    .label-group { opacity: 0; pointer-events: none; transition: opacity 0.1s ease-in-out; }')
    svg.append('    .det-group:hover .label-group { opacity: 1; }')
    
    svg.append(f'    .label-text {{ font-size: {font_size}px; font-weight: bold; fill: white; font-family: Arial, sans-serif; }}')
    svg.append('  </style>')
    
    for i, box in enumerate(boxes):
        x1, y1, x2, y2 = box["x1"], box["y1"], box["x2"], box["y2"]
        w_box = x2 - x1
        h_box = y2 - y1
        class_name = box["class_name"]
        conf = box["confidence"]
        c = box["class_id"]
        box_id = box["id"]                                  # <-- 取出 id

        b, g, r = colors(c, True)
        color_rgb = f"rgb({r},{g},{b})"  
        
        label_text = f"[{box_id}] {class_name} {conf:.1%}"  # <-- 标签前加 [id]
        
        svg.append(f'  <g class="det-group">')
        svg.append(f'    <rect class="det-box" x="{x1}" y="{y1}" width="{w_box}" height="{h_box}" pointer-events="all" />')
        svg.append('    <g class="label-group">')
        
        label_x = x1 + 2
        label_y = max(y1 - 5, font_size + 10)
        
        bg_width = len(label_text) * (font_size * 0.6) + font_size * 0.5
        bg_height = font_size + 8
        
        svg.append(f'      <rect x="{label_x}" y="{label_y - font_size}" width="{bg_width}" height="{bg_height}" rx="3" fill="{color_rgb}" />')
        svg.append(f'      <text class="label-text" x="{label_x + font_size * 0.2}" y="{label_y + 2}">{escape(label_text)}</text>')
        svg.append('    </g>')
        
        svg.append('  </g>')
    
    svg.append('</svg>')

    Path(svg_path).write_text("\n".join(svg), encoding="utf-8")
    print(f"SVG 已保存至: {svg_path}")


if __name__ == "__main__":
    import argparse
    # example usage: python predict.py --weights ./best.pt --source ./ANT_OR_14_01_62.JPG --conf 0.1 --output ./result.svg
    parser = argparse.ArgumentParser(description="YOLO Predict and Save Interactive SVG")
    parser.add_argument("--weights", type=str, default=r"./best.pt", help="Path to model weights")
    parser.add_argument("--source", type=str, default=r"./ANT_OR_14_01_62.JPG", help="Path to input image")
    parser.add_argument("--conf", type=float, default=0.1, help="Confidence threshold")
    parser.add_argument("--output", type=str, default="", help="Path to output SVG (default: auto-generated based on input name)")
    
    args = parser.parse_args()

    model_weights_path = args.weights
    input_img_dir = args.source
    prob_threshold = args.conf
    
    if args.output:
        output_svg_path = args.output
    else:
        # Default auto-generation if not explicitly provided
        output_svg_path = input_img_dir.split("/")[-1].split(".")[0] + f"_pred_{args.conf:.3f}.svg"


    model = YOLO(model_weights_path)

    results = model.predict(input_img_dir, conf=prob_threshold)

    for result in results:
        boxes = result.boxes 
        save_interactive_svg(result, output_svg_path)
    