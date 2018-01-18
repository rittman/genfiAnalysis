"""
Images for the GENFI paper
15th December 2016
"""

from os import path
from PIL import Image

# path definitions
out_directory = "figures"
wbvsaoo = "wholeBrainVsAOOResults"
wbvsaoo_lobar = "Reports/lobar_reports/wholeBrainVsAOOResults"
wbr_reports = "Reports/wholeBrainResults"
wbr = "wholeBrainResults"
wbr_lobes = "Reports/lobar_reports/wholeBrainResults"

# set image widths
#PNAS widths
c1w = 87 #mm
c2wa = 114 #mm
c2wb = 178 #mm

# set dpi
dpi = 600
ppm = float(dpi)/25.4 # convert to pixels per millimeter

def importResize(inList, outwd, ppmm):
    """
    This script takes two images and places them side by side with a fixed width
    (outw) measured in mm, with a pixels per millimeter values of ppm. The
    heights of the images are optimised so that both images fit the width.

    """
    wList = [] # list of image widths
    hList = [] # list of image heights
    imgList_in = [] # list of images

    # convert outwidth from mm to pixels using pixels per millimeter value
    outwd *= ppmm

    # ensure the out width is an integer
    outwd = int(outwd)

    for img in inList: # iterate through the input list
        im_in = Image.open(img) # import image
        w_in, h_in = im_in.size # get height and width
        wList.append(float(w_in)) # append width to list
        hList.append(float(h_in)) # append length to list
        imgList_in.append(im_in) # append image to list

    
    # output height
    max_h = max(hList)
    max_h_n = hList.index(max_h)
    
    for n,input_image in enumerate(imgList_in):
        wList[n]*=(max_h/hList[n])
        hList[n]*=(max_h/hList[n])
    # get ratio
    # multiplier for the first image
    total_width_input = float(sum(wList))

#    out_h = int((float(outwd)/float(wList[max_h_n])) * max_h)
    out_h = int(max_h * (outwd/(total_width_input/float(wList[max_h_n]))/wList[max_h_n]))
    print([outwd,out_h])
    
    outIm = Image.new("RGB", (outwd, out_h), "white") # create the output image
    
    left_a = 0
    upper_a = 0
    
    
    for n,input_image in enumerate(imgList_in):
        # percentage of total width
        w_pc = float(wList[n])/total_width_input
        
        # output width of panel
        outwd_panel = int(int(outwd) * w_pc)
        
        # output panel height
        outh_panel = int((float(outwd_panel)/float(wList[n]))*hList[n])
       
        
        
    
        # resize the image
        im_in = imgList_in[n].resize((outwd_panel,outh_panel))
        
        left_b = left_a
        upper_b = upper_a
        right_b = left_b+outwd_panel
        lower_b = upper_b+outh_panel
        
        outIm.paste(im_in, (left_b, upper_b, right_b, lower_b))
        
        left_a = right_b
        
    
#    y_r = ((hList[1]*outwd) / ((wList[0]*hList[1]) + (hList[0]*wList[1])))
#    
#    
#    x_r = ((y_r*hList[0]) / hList[1]) # multiplier for the second image
#
#    w0 = int(y_r*wList[0]) # new width for first image
#    h0 = int(y_r*hList[0]) # new height for first image
#
#    w1_a = int(x_r*wList[1]) # new width for second image
#    h1_a = int(x_r*hList[1]) # new height for second image
#
#    im0_in = imgList_in[0].resize((w0, h0)) # resize first image
#    im1_in = imgList_in[1].resize((w1_a, h1_a)) # resize second image
#
#    outIm = Image.new("RGB", (outwd, h0), "white") # create the output image
#
#    # paste the first image
#    left_a = 0
#    upper_a = 0
#    right_a = w0
#    lower_a = h0
#
#    outIm.paste(im0_in, (left_a, upper_a, right_a, lower_a))
#
#    # paste the first image
#    left_b = right_a
#    upper_b = 0
#    right_b = left_b+w1_a
#    lower_b = h1_a
#
#    outIm.paste(im1_in, (left_b, upper_b, right_b, lower_b))

    return outIm # return the output image

# Image 1 - top panel
img1 = path.join(wbr_reports, "degree_wt_allgroups.png")
#img2 = path.join(wbvsaoo, "degree_wtGenePos.png")
img2 = path.join(wbr_reports, "degree_wt_byGene.png")
 
img3 = path.join(wbvsaoo_lobar, "degree_wtGenePos_DiscontBkpt.png")

outImgTop = importResize([img1, img2, img3], c1w, ppm)

# Image 1 - bottom panel
img1 = path.join(wbr_reports, "geNorm_allgroups.png")
#img2 = path.join(wbvsaoo, "geNormGenePos.png")
img2 = path.join(wbr_reports, "geNorm_byGene.png")
img3 = path.join(wbvsaoo_lobar, "geNormGenePos_DiscontBkpt.png")

outImgBottom = importResize([img1, img2, img3], c1w, ppm)

# get the heights of the top and bottom panel
wt, ht = outImgTop.size
wb, hb = outImgBottom.size

# create the output image
outImg = Image.new("RGB", (wt, ht+hb))

# paste the top panel
left = 0
upper = 0
right = wt
lower = ht
outImg.paste(outImgTop, (left, upper, right, lower))

# paste the bottom panel
left = 0
upper = lower
right = wb
lower = lower + hb
outImg.paste(outImgBottom, (left, upper, right, lower))

# save the output as Figure1.png
outImg.info["dpi"] = dpi
outImg.save(path.join(out_directory, "Figure1.png"))

### Figure 4 - breakpoint analysis
img1 = path.join(wbvsaoo_lobar, "degree_wtGenePos_DiscontBkpt.png")
img2 = path.join(wbvsaoo_lobar, "geNormGenePos_DiscontBkpt.png")

im1 = Image.open(img1)
w1, h1 = im1.size

im2 = Image.open(img2)
w2, h2 = im2.size

# ratio of old to new lengths
rw = (c1w*ppm)/float(w1)

w1 = int(rw*w1)
h1 = int(rw*h1)
w2 = int(rw*w2)
h2 = int(rw*h2)

im1 = im1.resize((w1, h1))
im2 = im2.resize((w2, h2))

outw = int(c1w*ppm)
outImg = Image.new("RGB", (outw, h1+h2), "white")

left = 0
upper = 0
right = w1
lower = h1
outImg.paste(im1, (left, upper, right, lower))

left = 0
upper = lower
right = w2
lower += h2
outImg.paste(im2, (left, upper, right, lower))

outImg.info["dpi"] = dpi
outImg.save(path.join(out_directory, "Figure4.png"))

## Figure 2
# define list of lobes
lobeList = ["Frontal",
            "Temporal",
            "Parietal",
            "Occipital",
            "Cerebellum",
            "Hippocampus",
            "Cingulate",
            "Insula"#,
            # "Subcortical"
           ]

imggList = []
imgdbList = []
ncol = 4
outw = int(c2wa*ppm)

for l in lobeList:
    # open group comparison image
    imag = path.join(wbr_lobes, "degree_wt_"+l+"_allgroups.png")
    imgg = Image.open(imag)

    # open discontinous breakpoint analysis
    imdb = path.join(wbvsaoo, "degree_wtGenePos_DiscontBkpt_"+l+".png")
    imgdb = Image.open(imdb)

    w, h = imgg.size
    wr = (float(outw)/w)*(1/float(ncol))
    wo = int(w*wr)
    ho = int(h*wr)

    imgg = imgg.resize((wo, ho))

    w, h = imgdb.size
    wr = (float(outw)/w)*(1/float(ncol))
    wo = int(w*wr)
    ho = int(h*wr)

    imgdb = imgdb.resize((wo, ho))

    imggList.append(imgg)
    imgdbList.append(imgdb)

allLen = len(imggList)*2
hout = int((allLen/ncol + allLen%ncol) * ho)

# +50 to stop the bottom figures' axis labels being chopped off
outImg = Image.new("RGB", (outw, hout+50), "white")

count = 0
for n, imgg in enumerate(imggList):
    # group difference
    w, h = imgg.size
    left = int(float(w)*int(count%ncol))
    upper = int(float(h)*int(count/ncol))
    right = left+w
    lower = upper+h

    outImg.paste(imgg, (left, upper, right, lower))

    count += 1

    # discontinuous breakpoint analysis
    imgdb = imgdbList[n]
    w, h = imgdb.size
    left = int(float(w)*int(count%ncol))
    upper = int(float(h)*int(count/ncol))
    right = left+w
    lower = upper+h

    # build in slight adjustment to lower breakpoint analysis results
    adj = 25
    upper += adj
    lower += adj

    outImg.paste(imgdb, (left, upper, right, lower))
    count += 1


outImg.info["dpi"] = dpi
outImg.save(path.join(out_directory, "Figure2.png"))

## Figure 3
# define list of lobes
lobeList = ["Frontal",
            "Temporal",
            "Parietal",
            "Occipital",
            "Cerebellum",
            "Hippocampus",
            "Cingulate",
            "Insula"#,
            # "Subcortical
           ]

imggList = []
imgdbList = []
ncol = 4
outw = int(c2wa*ppm)

for l in lobeList:
    # open group comparison image
    imag = path.join(wbr_lobes, "closeCentNorm_"+l+"_allgroups.png")
    imgg = Image.open(imag)

    # open discontinous breakpoint analysis
    imdb = path.join(wbvsaoo, "closeCentNormGenePos_DiscontBkpt_"+l+".png")
    imgdb = Image.open(imdb)

    w, h = imgg.size
    wr = (float(outw)/w)*(1/float(ncol))
    wo = int(w*wr)
    ho = int(h*wr)

    imgg = imgg.resize((wo, ho))

    w, h = imgdb.size
    wr = (float(outw)/w)*(1/float(ncol))
    wo = int(w*wr)
    ho = int(h*wr)

    imgdb = imgdb.resize((wo, ho))

    imggList.append(imgg)
    imgdbList.append(imgdb)

allLen = len(imggList)*2
hout = int((allLen/ncol + allLen%ncol) * ho)
outImg = Image.new("RGB", (outw, hout), "white")

count = 0
for n, imgg in enumerate(imggList):
    w, h = imgg.size
    left = int(float(w)*int(count%ncol))
    upper = int(float(h)*int(count/ncol))
    right = left+w
    lower = upper+h

    outImg.paste(imgg, (left, upper, right, lower))

    count += 1

    imgdb = imgdbList[n]
    w, h = imgdb.size
    left = int(float(w)*int(count%ncol))
    upper = int(float(h)*int(count/ncol))
    right = left+w
    lower = upper+h

    # build in slight adjustment to lower breakpoint analysis results
    adj = 25
    upper += adj
    lower += adj

    outImg.paste(imgdb, (left, upper, right, lower))
    count += 1


outImg.info["dpi"] = dpi
outImg.save(path.join(out_directory, "Figure3.png"))


# Figure 5 - hubs
img1 = path.join(wbr, "degree_wt_allgroups.png")
img2 = path.join(wbvsaoo, "degree_wtGenePos.png")

outImgTop = importResize([img1, img2], c1w, ppm)

# Figure 6 - volume vs connectivity
imgList = ["volVsConnectivity/volumesVsConnectivity_Whole Brain Volume.png",
           "volVsConnectivity/volumesVsConnectivity_Frontal lobe volume.png",
           "volVsConnectivity/volumesVsConnectivity_Temporal lobe Volume.png",
           "volVsConnectivity/volumesVsConnectivity_Parietal lobe Volume.png"] #,
#           "volVsConnectivity/volumesVsGlobalEfficiency_Whole Brain Volume.png",
#           "volVsConnectivity/volumesVsGlobalEfficiency_Frontal lobe volume.png",
#           "volVsConnectivity/volumesVsGlobalEfficiency_Temporal lobe Volume.png",
#           "volVsConnectivity/volumesVsGlobalEfficiency_Parietal lobe Volume.png"]

outImg = None

nrow=2
ncol = len(imgList) / nrow

for n, i in enumerate(imgList):
    im = Image.open(i)
    if not outImg:
        w, h = im.size
        outw = int(c1w*ppm)
    
        # multiplier for the images - assumes images are the same size
        y = ((float(h)*outw) / ((float(w)*h)*ncol))

    # resize images
    imw = int(y*w)
    imh = int(y*h)

    im = im.resize((imw, imh))

    if not outImg:
        outImg = Image.new("RGB", (outw, imh*nrow), "white")

    # paste on images
    left = int(n%ncol * imw)
    upper = int(n/ncol) * imh
    right = left+imw
    lower = upper+imh

    outImg.paste(im, (left, upper, right, lower))

outImg.save(path.join(out_directory, "Figure5.png"))
