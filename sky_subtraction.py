import numpy as np
from astropy.io import fits
import pylab as plt
import glob as glob
from matplotlib.colors import LogNorm
from scipy.stats import mode
from specutils import get_wcs_solution
from scipy.interpolate import interp1d
from mask import mask
from center_spec_row import finding_row_center

#data1 = fits.getdata('mods1b.20161118.0017_otfc_derot.fits.gz')
#data2 = fits.getdata('mods2b.20161118.0015_otfc_derot.fits.gz')

#row_center_guess = 1690

#cut_data1_below = data1[row_center - 160: row_center - 40,:]
#cut_data1_above = data1[row_center + 40: row_center + 160,:]


#cut_data2_below = data2[row_center - 160: row_center - 40,:]
#cut_data2_above = data2[row_center + 40: row_center + 160,:]

#split_data1_below = np.vsplit(cut_data1_below, 4)
#split_data1_above = np.vsplit(cut_data1_above, 4)


#split_data2_below = np.vsplit(cut_data2_below, 4)
#split_data2_above = np.vsplit(cut_data2_above, 4)

files = [x for x in glob.glob('*.gz')]

def myround(x, base=5):
    return int(base * round(x/base))


def sky_sum(data, wavelength, spec_wavelength):
    
    '''
    What this function does is that you pass in a chunk of data from above or below your taget spectrum plus or minus a window 
    so as to not include your spectrum and this function will find the sky spectrum from the data you pass in.

    Parameter
    -------------
    data: this is the chunk of 2D spectrum data that will be analyzed to calculate the sky spectrum. 
          Make sure that the input data does not have any
          spectrum data in it as this will throw off your sky subtraction and the length of the row is divisible by 5.
    
    wavelength: this is a 2D array in the same portion as your data array and we need this to implement interpolation 

    spec_wavelength: This is a 1D array/list that holds the values of the wavelength of the target spectrum.
                     This is used to map the values from interpolation onto the wavelength of our target spectrum

    Outputs
    -------------
    sky_spec: a 1D array of spectrum correspoding to the data passed in.
    
    '''
    

    #splitting the data chunk into 5 pieces row-wise
    split_data = np.vsplit(data, 5)
    
    #splitting the wavelength array into 5-pieces row-wise
    split_wave = np.vsplit(wavelength, 5)
    
    #making a list that will hold the arrays of each chunk of sky, at the end this should be len=5 as im
    #splitting the array into 5 chunks
    sky_chunks = []
    
    #for loop that loops through each row and finds the intepoation and later sums up the sky spectrum from 
    for i in range(len(split_data)):
        
        #variable that will hold the sum of the sky for each of the 5 chunks
        sum_sky = 0
        
        #variable to keep track how many times we looped through the row and then divide sum_sky by N to get the average
        N = 0

        for j in range(len(split_data[i][:,0])):
            
            #getting the interpolation for our sky
            f = interp1d(split_wave[i][j,:], split_data[i][j,:], fill_value='extrapolate')

            #getting the sum_sky spectrum by taking the interpolation and assing in the spectrum wavelength as arguments
            #we do this to map the sky onto our target spectrum
            sum_sky += f(spec_wavelength)
            
            #incrementing counter by 1
            N+=1
        
        sky_chunks.append(sum_sky/N)

    return sky_chunks

def test_sky(sky_chunk, above): 
    
    '''
    This function will check to see if the sky changes the farther away we go from the target spectrum.

    Parameter
    --------------
    sky_chunk: this is the output from the function sky_sum and it should be chunks of sky away from the target
    above: this is a boolean variable and will dicitate which test to take in the function below
    
    Output:
    med_sky: A median sky spectrum of good sky near the target
    
    '''
    
    if above == True:
        
        #this takes the difference between the spectrums, this is to check if the sky varies as we move away from the source 
        diff = sky_chunk - sky_chunk[0]
        
        #making list that will hold what is the number that occurs the most in the diff arrays 
        num_most = []
    
        for i in diff:
        
            most = mode(i, axis=None)[0][0]
            
            #print(most)   
            
            if most < 2:   
                num_most.append(True)
            
            else:
                num_most.append(False)
        
        #making a sky spectrum
        med_sky = np.median(np.array(sky_chunk)[num_most], axis = 0)
        
        #print(len(num_most))

        return med_sky
        
    if above == False:
        
        #this takes the difference between the spectrums, this is to check if the sky varies as we move away from the source 
        diff = sky_chunk - sky_chunk[-1]
        
        
        #making list that will hold what is the number that occurs the most in the diff arrays 
        num_most = []
        
        for i in diff:
        
            most = mode(i, axis=None)[0][0]
            
            #print(most)

            if most < 2:   
                num_most.append(True)
            
            else:
                num_most.append(False)

        #making median sky spectrum
        #print(len(num_most))
        med_sky = np.median(np.array(sky_chunk)[num_most], axis = 0)

        return med_sky

def master_sky(sky_above, sky_below):
    
    '''
    This function will return the master sky spectrum by taking the median of sky_above and sky_below

    Parameters
    -------------------
    sky_above: this is the 1D sky spec above the target and is the output from sky_sum using the above data
    sky_below: this is the 1D sky spec below the target and is the output from sky_sum using the below data

    Output
    -------------------
    sky_master: this is the master sky spectrum we will use to subtract the sky form every row in the data array.

    '''
    sky_master =  np.median([sky_above, sky_below], axis = 0)
    
    return sky_master
    
    
    
    
    
def sky_spec(filename):
    
    '''
    This function will try to find the spectrum of the sky given a liist/array
    that holds the 2D information of the sky. Note that this is already split meaning
    that we have a list/aray with each element representing all the columns and 
    certain rows
    
    Parameter
    -------------
    data: This is the array/list of data to get the sky spectrum

    Output
    -------------
    A Median sky spectrum

    '''
    
    #getting the data from the filename
    data = fits.getdata(filename)
    
    #print('filename: ' + filename+ ', shape =  '+ str(data.shape))

    #getting the header for the file as we will need it to get transformation from pixels to wavelength.
    hdr = fits.getheader(filename)
    
    #making an array with pixel-coordinates
    pix_x, pix_y = np.meshgrid(range(len(data[0,:])), range(len(data[:,0])))
    
    #getting the transormation for the file from the header
    transformation = get_wcs_solution(hdr)
    
    #this gives me a 2D wavelength array for each (x, y) pixel coordinates in the array
    wvln_array = transformation(pix_x, pix_y)

    #making my guess to the row center into its own variable
    row_center_guess = 1690
    
    #making my max and minimum variables which will limit the box around what I want 
    minimum = 200
    maximum = 125

    #making the actual row_min variables and row_max variables
    row_min = row_center_guess - minimum
    row_max = row_center_guess + maximum
    
    #getting the mask that will filter out the bad data
    filt_mask = mask(row_min, row_max, data)
    
    #applying the mask to the data
    reduced_data = filt_mask * data
    
    #getting the index of the row for the center of our spectrum
    spec_row_center = finding_row_center(reduced_data)

    #rounding the index to the nearest 5 or 0 so that vsplit can work properly in the sum_sky function
    row_center = myround(spec_row_center)
    
    #applying mask to wavelength array
    wvln_reduced = wvln_array * filt_mask

    #making a window to not include the target spectrum in this case we go 50 pixels above and below
    window = 50

    #making variables that will hold data above and below the window so as to find the sky spectrum
    data_below = reduced_data[row_min: row_center - window,:]
    data_above = reduced_data[row_center + window: row_max :]

    #similar for the wavelength as we will need this for sum_sky function
    wvln_below = wvln_reduced[row_min: row_center - window,:]
    wvln_above = wvln_reduced[row_center + window: row_max :]

    #getting the wavelength information of where the center row of our spectrum is from the function above
    spec_wvln = wvln_reduced[spec_row_center]
    
    sky_chunk_above = sky_sum(data_above, wvln_above, spec_wvln)
    sky_chunk_below = sky_sum(data_below, wvln_below, spec_wvln)
    
    sky_above = test_sky(sky_chunk_above, True)
    sky_below = test_sky(sky_chunk_below, False)


    sky_master = master_sky(sky_above, sky_below)
    
    test = data - sky_master
    
    #plt.figure(figsize = (10,10))
    #plt.imshow(test, origin='lower', cmap ='gray', norm=LogNorm())
    #plt.show()
    
    return test

avg_spec = 0
m = 0
for i in files:
    
    sky_spec(i)
    m+=1

avg = avg_spec/m

plt.figure(figsize = (10,10))
plt.imshow(avg, origin='lower', cmap ='gray', norm=LogNorm())
plt.show()
#sky_spec(split_data1_below)
#sky_spec(split_data1_above)
#sky_spec(split_data2_below)
#sky_spec(split_data2_above)


