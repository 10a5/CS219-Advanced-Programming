import cv2
import argparse
import sys
import os
import numpy as np

def adjust_brightness(img, delta):
    lut = np.clip(np.arange(256, dtype=np.int16) + delta, 0, 255).astype(np.uint8)
    return cv2.LUT(img, lut)

def blend_images(img1, img2):
    if img1.shape[:2] != img2.shape[:2]:
        img2 = cv2.resize(img2, (img1.shape[1], img1.shape[0]))
    if img1.ndim != img2.ndim:
        if img1.ndim == 2:
            img1 = cv2.cvtColor(img1, cv2.COLOR_GRAY2BGR)
        else:
            img2 = cv2.cvtColor(img2, cv2.COLOR_GRAY2BGR)
    
    return cv2.addWeighted(img1, 0.5, img2, 0.5, 0)

def validate_input(args):
    valid_ops = ["brightness", "blend"]
    if args.operation not in valid_ops:
        print(f"Error: Invalid {args.operation}, Aviliable: {valid_ops}")
        sys.exit(1)
    
    required_images = 1 if args.operation == "brightness" else 2
    if len(args.images) != required_images:
        print(f"Error: {args.operation} operation requires {required_images} images, detected {len(args.images)}")
        sys.exit(1)
    
    for path in args.images:
        if not os.path.isfile(path):
            print(f"Error: Path doesn't exist [{path}]")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Image Processing Tool")
    parser.add_argument("-i", "--images", nargs='+', required=True,
                      help="Enter image paths (space-separated)")
    parser.add_argument("-op", "--operation", required=True,
                      choices=["brightness", "blend"],
                      help="Operation type: brightness or blend")
    parser.add_argument("-b", "--brightness", type=int, default=30,
                      help="Brightness adjustment value (default 30)")
    parser.add_argument("-o", "--output", default="output.jpg",
                      help="Output image path (default output.jpg)")
    
    args = parser.parse_args()
    validate_input(args)
    images = [cv2.imread(path) for path in args.images]
    if args.operation == "brightness":
        result = adjust_brightness(images[0], args.brightness)
    else:
        result = blend_images(*images)

    cv2.imshow("Result", result)
    print("Press to exit...")
    cv2.waitKey(0) 
    cv2.destroyAllWindows()

    cv2.imwrite(args.output, result)
    print(f"Result save to: {args.output}")

if __name__ == "__main__":
    main()