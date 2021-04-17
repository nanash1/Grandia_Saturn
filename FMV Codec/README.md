# Saturn Grandia FMV codec

Reverse engineered Saturn Grandia FMV codec. Written in Python 3 with Numpy. This is WIP.

## Usage

The usage is straight forward with some working knowledge of Python. Install the package by changing to the directory that contains setup.py and run:
```
pip install .
```
The following example code shows how to use the package. You could, of course, use a different method to generate the required numpy array, but I found this to be the fastest.

```python
import numpy as np
import ss_grandia_codec
from decord import VideoReader
from decord import cpu
from PIL import Image

compression_level = 10
frame_number = 900

# Import video with decord's VideoReader
frame = VideoReader("examples/intro.mkv", ctx=cpu(0))[frame_number].asnumpy()

# Read into PIL.Image for downscaling and croping
frame = Image.fromarray(frame, mode='RGB')
frame = frame.resize((352, 198), resample=Image.LANCZOS)
frame = frame.crop((0, 11, 352, 187))

# Convert to numpy RGB array
frame = np.array(frame)

# Encode frame to raw video data
frame_data = ss_grandia_codec.encode_frame(frame, compression_level, subsample="avrg", yscale=(16,235), cscale=224)

# Write the result to a file
with open("encoded.bin", "wb") as outf:
    outf.write(frame_data)
```

The encode_frame function has some optional keyword arguments, so make sure to read the docstring.
