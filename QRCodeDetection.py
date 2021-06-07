
from matplotlib import pyplot
from matplotlib.patches import Rectangle
from pyzbar.pyzbar import decode

import numpy as np
import imageIO.png
import math

class Queue:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)


def createInitializedGreyscalePixelArray(image_width, image_height, initValue = 0):

    new_array = [[initValue for x in range(image_width)] for y in range(image_height)]
    return new_array


# this function reads an RGB color png file and returns width, height, as well as pixel arrays for r,g,b
def readRGBImageToSeparatePixelArrays(input_filename):

    image_reader = imageIO.png.Reader(filename=input_filename)
    # png reader gives us width and height, as well as RGB data in image_rows (a list of rows of RGB triplets)
    (image_width, image_height, rgb_image_rows, rgb_image_info) = image_reader.read()

    print("read image width={}, height={}".format(image_width, image_height))

    # our pixel arrays are lists of lists, where each inner list stores one row of greyscale pixels
    pixel_array_r = []
    pixel_array_g = []
    pixel_array_b = []

    for row in rgb_image_rows:
        pixel_row_r = []
        pixel_row_g = []
        pixel_row_b = []
        r = 0
        g = 0
        b = 0
        for elem in range(len(row)):
            # RGB triplets are stored consecutively in image_rows
            if elem % 3 == 0:
                r = row[elem]
            elif elem % 3 == 1:
                g = row[elem]
            else:
                b = row[elem]
                pixel_row_r.append(r)
                pixel_row_g.append(g)
                pixel_row_b.append(b)

        pixel_array_r.append(pixel_row_r)
        pixel_array_g.append(pixel_row_g)
        pixel_array_b.append(pixel_row_b)

    return (image_width, image_height, pixel_array_r, pixel_array_g, pixel_array_b)

# This method packs together three individual pixel arrays for r, g and b values into a single array that is fit for
# use in matplotlib's imshow method
def prepareRGBImageForImshowFromIndividualArrays(r,g,b,w,h):
    rgbImage = []
    for y in range(h):
        row = []
        for x in range(w):
            triple = []
            triple.append(r[y][x])
            triple.append(g[y][x])
            triple.append(b[y][x])
            row.append(triple)
        rgbImage.append(row)
    return rgbImage
    

# This method takes a greyscale pixel array and writes it into a png file
def writeGreyscalePixelArraytoPNG(output_filename, pixel_array, image_width, image_height):
    # now write the pixel array as a greyscale png
    file = open(output_filename, 'wb')  # binary mode is important
    writer = imageIO.png.Writer(image_width, image_height, greyscale=True)
    writer.write(file, pixel_array)
    file.close()

#rgb to greyscale
def computeRGBToGreyscale(pixel_array_r, pixel_array_g, pixel_array_b, image_width, image_height):  
    greyscale_pixel_array = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range (image_height):
        for j in range(image_width):
            greyscale_pixel_array[i][j] = round(.299*pixel_array_r[i][j] + .587*pixel_array_g[i][j] + .114*pixel_array_b[i][j])    
    
    
    return greyscale_pixel_array
#scale array between 0-255
def scaleTo0And255AndQuantize(pixel_array, image_width, image_height):
    max = 0
    min = 99999
    for i in range(len(pixel_array)):
        for value in range(len(pixel_array[i])):
            pixel = pixel_array[i][value]
            if pixel > max:
                max = pixel
            if min > pixel:
                min = pixel
    
    r = (max-min)
    if r == 0:
        c = 0
    else:
        c = 255/(max-min)
    for i in range(len(pixel_array)):
        for value in range(len(pixel_array[i])):
            p = pixel_array[i][value]
            pixel_array[i][value] = round(c*(p-min))
    return pixel_array

def computeHorizontalEdgesSobel(pixel_array, image_width, image_height):
    r = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(1, image_height-1):
        for j in range(1, image_width-1):
            r[i][j] = (pixel_array[i-1][j-1] + (pixel_array[i-1][j] * 2) + pixel_array[i-1][j+1] + (pixel_array[i+1][j-1] * -1) + (pixel_array[i+1][j] * -2) + (pixel_array[i+1][j+1] * -1))/8
    return r

def computeVerticalEdgesSobel(pixel_array, image_width, image_height):
    greyscale_pixel_array = createInitializedGreyscalePixelArray(image_width, image_height)
    pixel_scale=[-1,0,1,-2,0,2,-1,0,1]
    for i in range(0,image_height-2):
        for j in range(0,image_width-2):
            greyscale_pixel_array[i+1][j+1] = pixel_array[i][j]*pixel_scale[0]+pixel_array[i][j+1]*pixel_scale[1]+pixel_array[i][j+2]*pixel_scale[2]+pixel_array[i+1][j]*pixel_scale[3]+pixel_array[i+1][j+1]*pixel_scale[4]+pixel_array[i+1][j+2]*pixel_scale[5]+pixel_array[i+2][j]*pixel_scale[6]+pixel_array[i+2][j+1]*pixel_scale[7]+pixel_array[i+2][j+2]*pixel_scale[8]
            greyscale_pixel_array[i+1][j+1] = (greyscale_pixel_array[i+1][j+1]/8)
    return greyscale_pixel_array

def computeGradiantMagnitude(vertical, horizontal, image_width, image_height):
    magnitude_array = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(1,image_height-1):
        for j in range(1, image_width-1):
            magnitude_array[i][j] = math.sqrt((vertical[i][j]*vertical[i][j]) + (horizontal[i][j]*horizontal[i][j]))
    return magnitude_array

def computeGaussianAveraging3x3RepeatBorder(pixel_array, image_width, image_height):
    gaussian_array = createInitializedGreyscalePixelArray(image_width, image_height)
    for a in pixel_array:
        a.insert(0, a[0])
        a.insert(-1, a[-1])
    pixel_array.insert(0, pixel_array[0])
    pixel_array.insert(-1, pixel_array[-1])
    for i in range(1, image_height+1):
        for j in range(1, image_width+1):
            gaussian_array[i-1][j-1] = (pixel_array[i-1][j-1] + (pixel_array[i-1][j]*2) + pixel_array[i-1][j+1] + (pixel_array[i][j-1]*2) + (pixel_array[i][j]*4)+ (pixel_array[i][j+1]*2) + pixel_array[i+1][j-1] + (pixel_array[i+1][j]*2) + pixel_array[i+1][j+1])/16    
    return gaussian_array

def computeThresholdGE(pixel_array, threshold_value, image_width, image_height):
    threshold_array = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            if pixel_array[i][j] < threshold_value:
                threshold_array[i][j] = 0
            else:
                threshold_array[i][j] = 255
    return threshold_array

def computeDilation8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    dilation = createInitializedGreyscalePixelArray(image_width, image_height)
    for array in pixel_array:
        array.insert(0, 0)
        array.append(0)
    zeros = [0] * (image_width + 2)
    pixel_array.insert(0, zeros)
    pixel_array.append(zeros)

    for i in range(1, image_height+1):
        for j in range(1, image_width+1):
            value = pixel_array[i-1][j-1] + pixel_array[i-1][j] + pixel_array[i-1][j+1] + pixel_array[i][j-1] + pixel_array[i][j] + pixel_array[i][j+1] + pixel_array[i+1][j-1] + pixel_array[i+1][j] + pixel_array[i+1][j+1]
            if value > 0:
                dilation[i-1][j-1] = 1
            else:
                dilation[i-1][j-1] = 0
    return dilation


def computeErosion8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    Erosion = createInitializedGreyscalePixelArray(image_width, image_height)
    for array in pixel_array:
        array.insert(0, 0)
        array.append(0)
    zeros = [0] * (image_width + 2)
    pixel_array.insert(0, zeros)
    pixel_array.append(zeros)

    for i in range(1, image_height+1):
        for j in range(1, image_width+1):
            if (pixel_array[i-1][j-1] != 0 and pixel_array[i-1][j] != 0 and pixel_array[i-1][j+1] != 0 and pixel_array[i][j-1] != 0 and pixel_array[i][j] != 0 and pixel_array[i][j+1] != 0 and pixel_array[i+1][j-1] != 0 and pixel_array[i+1][j] != 0 and pixel_array[i+1][j+1] != 0):
                Erosion[i-1][j-1] = 1
            else:
                Erosion[i-1][j-1] = 0

    return Erosion


def computeConnectedComponentLabeling(pixel_array, image_width, image_height):
    connected = createInitializedGreyscalePixelArray(image_width, image_height)
    visited = createInitializedGreyscalePixelArray(image_width, image_height)
    setter = 1
    dictionary = dict()
    Q = Queue()
    
    for i in range(1, image_height):
        for j in range(1, image_width):
            if pixel_array[i][j] and not visited[i][j]:
                
                Q.enqueue([i,j])
                visited[i][j] = True
                
                while not Q.isEmpty():
                    head = Q.dequeue()
                    if setter not in dictionary:
                        dictionary[setter] = 1
                        connected[head[0]][head[1]] = setter
                    else:
                        dictionary[setter] += 1
                        connected[head[0]][head[1]] = setter
                            
                    if head[0]-1 >= 0 and not visited[head[0]-1][head[1]] and pixel_array[head[0]-1][head[1]]:
                        Q.enqueue([head[0]-1,head[1]])
                        visited[head[0]-1][head[1]] = True
                        
                    if  head[1]-1 >= 0 and not visited[head[0]][head[1]-1] and pixel_array[head[0]][head[1]-1]:
                        Q.enqueue([head[0],head[1]-1])
                        visited[head[0]][head[1]-1] = True
                        
                    if head[0]+1 < image_height and not visited[head[0]+1][head[1]] and pixel_array[head[0]+1][head[1]]:
                        Q.enqueue([head[0]+1,head[1]])
                        visited[head[0]+1][head[1]] = True
                        
                    if  head[1]+1 < image_width and not visited[head[0]][head[1]+1] and pixel_array[head[0]][head[1]+1]:
                        Q.enqueue([head[0],head[1]+1])
                        visited[head[0]][head[1]+1] = True
                setter += 1  
                
    return connected, dictionary

def CleanArray(pixel_array, image_width, image_height, max_key):
    clean = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            if pixel_array[i][j] != max_key:
                clean[i][j] = 0
            else:
                clean[i][j] = 255
    return clean

def FindCorner(pixel_array, image_width, image_height):
    min_x = 99999
    max_x = -10
    min_y = 99999
    max_y = -10
    for i in range(image_height):
        for j in range(image_width):
            if pixel_array[i][j] > 0:
                if j > max_x:
                    max_x = j
                if j < min_x:
                    min_x = j
                if i > max_y:
                    max_y = i    
                if i < min_y:
                    min_y = i
    return min_x, max_x, min_y, max_y          

def GetQRpixel(pixel_array, min_x, max_x, min_y, max_y):
    qr = []
    width = 0
    height = 0
    for i in range(min_y, max_y):
        row = []
        height += 1
        width = 0
        for j in range(min_x, max_x):
            width += 1
            row.append(pixel_array[i][j])
        qr.append(row)
    return qr, width, height

def main():
    filename = "./images/covid19QRCode/poster1small.png"
    #filename = "./images/covid19QRCode/challenging/bch.png"
    #filename = "./images/covid19QRCode/challenging/bloomfield.png"


    # we read in the png file, and receive three pixel arrays for red, green and blue components, respectively
    # each pixel array contains 8 bit integer values between 0 and 255 encoding the color values
    (image_width, image_height, px_array_r, px_array_g, px_array_b) = readRGBImageToSeparatePixelArrays(filename)
    original = prepareRGBImageForImshowFromIndividualArrays(px_array_r, px_array_g, px_array_b, image_width, image_height)

    #convert rgb to greyscale and scale it
    greyscale = computeRGBToGreyscale(px_array_r, px_array_g, px_array_b,image_width, image_height)
    greyscale_scaled = scaleTo0And255AndQuantize(greyscale,image_width,image_height)

    #horizontal and vertical edge
    horizontal = computeHorizontalEdgesSobel(greyscale_scaled,image_width,image_height)
    vertical = computeVerticalEdgesSobel(greyscale_scaled,image_width,image_height)

    #compute edge gradiant for the image
    mag = computeGradiantMagnitude(vertical,horizontal,image_width, image_height)
    
    #Smooth over the edge mag with gaussian
    gaussian = mag
    for i in range(8):
        gaussian = computeGaussianAveraging3x3RepeatBorder(gaussian,image_width,image_height)
    gaussian_stretch = scaleTo0And255AndQuantize(gaussian,image_width,image_height)


    #Threshold
    thresh = computeThresholdGE(gaussian_stretch, 80, image_width, image_height)

    #findhole
    dilation = thresh
    for i in range(2):
        dilation = computeDilation8Nbh3x3FlatSE(dilation,image_width,image_height)

    #connect components
    print("finding connected")
    connected, size = computeConnectedComponentLabeling(dilation,image_width,image_height)
    max_key = max(size, key=size.get)
    clean = CleanArray(connected, image_width,image_height, max_key)

    #Finding Coordinate
    min_x, max_x, min_y, max_y = FindCorner(clean,image_width,image_height)
    rect_width = max_x - min_x
    rect_height = max_y - min_y
    
    #Getting Just QR Code (EXTENSION)
    QR, width, height = GetQRpixel(original, min_x, max_x, min_y, max_y)
    array = np.array(QR)
    pyplot.imsave('QR.png', array)
    #setplot figure
    pyplot.imshow(original, cmap='gray')

    # get access to the current pyplot figure
    axes = pyplot.gca()
    # create a 70x50 rectangle that starts at location 10,30, with a line width of 3
    rect = Rectangle( (min_x, min_y), rect_width, rect_height, linewidth=3, edgecolor='g', facecolor='none' )
    # paint the rectangle over the current plot
    axes.add_patch(rect)

    # plot the current figure
    pyplot.show()



if __name__ == "__main__":
    main()