import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.image as mimg

img = mimg.imread('../../data/png/distr_bf97.png')

fig, ax = plt.subplots()
rect = ptc.Rectangle((0.25, 0.3), width=0.5, height=0.4, fill=False)
ax.add_patch(rect)
ax.imshow(img, extent=(0.25, 0.75, 0.3, 0.7))
ax.axis('off')
ax.set_xlim((0,1))
ax.set_ylim((0,1))
ax.text(0.5, 0.25, 'Модель', ha='center', va='center')
plt.show()