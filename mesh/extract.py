#%%
import numpy as np
from PIL import Image, ImageFilter
from matplotlib import pyplot as plt

img = Image.open('skull.png').convert('L')
img = Image.open('tree.png').convert('L')
img = Image.open('circle.jpg').convert('L')
img = Image.open('rose.webp').convert('L')
light = np.array(img.filter(ImageFilter.GaussianBlur(1.5)))
light = light[::-1]/light.max()

# Real distance between pixels [m] (replace with your own value)
dx = 1/np.sqrt(light.size)

contour = plt.contour(light, levels=[0.5], colors='r')
paths = contour.collections[0].get_paths()
plt.close()
# %%
curves = []
for path in paths:
    vertecies = path.vertices 
    vertecies = vertecies.round(4) * dx
    start = vertecies[0] 
    
    i = 1
    while(i < len(vertecies)):
        if np.all(vertecies[i] == start):
            curves.append(vertecies[:i+1])
            if i+1 < len(vertecies):
                vertecies = vertecies[i+1:]
                start = vertecies[0]
                i = 1
            else:
                break
        else:
            i += 1

print(f'found {len(curves)} curves')
for c in curves:
    plt.plot(c[:,0], c[:,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#%%

np.save('curves.npy', np.array(curves, dtype=object))

# %%
data = np.load('curves.npy', allow_pickle=True)
