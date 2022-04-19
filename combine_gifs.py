import imageio
import numpy as np
import os
import numpy as np
import PIL
from PIL import Image
import re

gifs = []

location = "2022_04_19_13_31_16_ARC_EUCLID_SQRT"
path_to_experiments = r'C:\\Users\\20175326\\Desktop\\Thesis\\Code\\height-based-flow-maps-java\\Height based flow maps\\experiments\\'  

valid_images = [".gif"]
arr = os.listdir(path_to_experiments + location)

imgs = []
gif_location = []
for item in arr:
    path = path_to_experiments + location
    path = path + "\\" + item
    for f in os.listdir(path):
        ext = os.path.splitext(f)[1]
        if ext.lower() not in valid_images:
            continue
        gif_location.append(os.path.join(path,f))

#gif_location   

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]
gif_location.sort(key=natural_keys)

gifs = []

for item in gif_location:
    gifs.append(imageio.get_reader(item))

number_of_frames = gifs[0].get_length()

new_gif = imageio.get_writer(location + '.gif', fps=1)

for frame_number in range(number_of_frames):
    counter = 0
    
    img1 = gifs[0].get_next_data()
    img2 = gifs[1].get_next_data()    
   
    new_image = np.hstack((img1, img2))
    
    for gif in gifs[2:]:
    
        img = gif.get_next_data()
        new_image = np.hstack((new_image, img))
    new_gif.append_data(new_image)


for gif in gifs:
    gif.close()
new_gif.close()

print("done!")