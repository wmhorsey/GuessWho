#!/usr/bin/env python3
"""
Render PLY frame sequence to PNG images using matplotlib.
Usage: python render_frames.py [frames_dir] [output_dir]

Then convert to video with:
  ffmpeg -framerate 30 -i renders/frame_%05d.png -c:v libx264 -pix_fmt yuv420p output.mp4
"""

import os
import sys
import struct

def read_ply(filename):
    """Read a simple ASCII PLY point cloud with RGB colors."""
    points = []
    colors = []
    
    with open(filename, 'r') as f:
        # Skip header
        line = f.readline()
        vertex_count = 0
        while line.strip() != 'end_header':
            if line.startswith('element vertex'):
                vertex_count = int(line.split()[-1])
            line = f.readline()
        
        # Read vertices
        for _ in range(vertex_count):
            parts = f.readline().split()
            if len(parts) >= 6:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                r, g, b = int(parts[3]), int(parts[4]), int(parts[5])
                points.append((x, y, z))
                colors.append((r/255.0, g/255.0, b/255.0))
    
    return points, colors

def render_frame(ply_file, output_file, azim=45, elev=20):
    """Render a PLY file to a PNG image."""
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    points, colors = read_ply(ply_file)
    
    if not points:
        print(f"  Warning: No points in {ply_file}")
        return False
    
    xs, ys, zs = zip(*points)
    
    fig = plt.figure(figsize=(10, 10), facecolor='black')
    ax = fig.add_subplot(111, projection='3d', facecolor='black')
    
    # Plot points
    ax.scatter(xs, ys, zs, c=colors, s=20, alpha=0.8)
    
    # Set viewing angle
    ax.view_init(elev=elev, azim=azim)
    
    # Style
    ax.set_xlim(min(xs)-10, max(xs)+10)
    ax.set_ylim(min(ys)-10, max(ys)+10)
    ax.set_zlim(min(zs)-10, max(zs)+10)
    ax.set_axis_off()
    
    # Add frame number to title
    frame_num = os.path.basename(ply_file).replace('frame_', '').replace('.ply', '')
    ax.set_title(f'Frame {frame_num}', color='white', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, facecolor='black', edgecolor='none')
    plt.close()
    
    return True

def main():
    frames_dir = sys.argv[1] if len(sys.argv) > 1 else 'frames'
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'renders'
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all PLY files sorted
    ply_files = sorted([f for f in os.listdir(frames_dir) if f.endswith('.ply')])
    
    if not ply_files:
        print(f"No PLY files found in {frames_dir}/")
        return
    
    print(f"Rendering {len(ply_files)} frames from {frames_dir}/ to {output_dir}/")
    
    for i, ply_file in enumerate(ply_files):
        input_path = os.path.join(frames_dir, ply_file)
        output_path = os.path.join(output_dir, ply_file.replace('.ply', '.png'))
        
        # Slowly rotate camera as animation progresses
        azim = 45 + (i * 2)  # Rotate 2 degrees per frame
        
        render_frame(input_path, output_path, azim=azim)
        
        if (i + 1) % 10 == 0 or i == len(ply_files) - 1:
            print(f"  Rendered {i + 1}/{len(ply_files)} frames")
    
    print(f"\n✅ Done! Images saved to {output_dir}/")
    print(f"\n💡 To create video, run:")
    print(f"   ffmpeg -framerate 30 -i {output_dir}/frame_%05d.png -c:v libx264 -pix_fmt yuv420p vortex_evolution.mp4")

if __name__ == '__main__':
    main()
