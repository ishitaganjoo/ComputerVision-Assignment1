//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>

using namespace std;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input);

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);

SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());

  cout<<"In the convolve general\n";
  int HEIGHT = input.rows();
  int WIDTH = input.cols();
  for(int row = 0; row < HEIGHT; row++)
  { 
  	for(int col = 0; col < WIDTH; col++ ) 
	{
  		float accumulation = 0;
  		float weightsum = 0; 
  		for(int i = -1; i <= 1; i++ )
		{
    			for(int j = -1; j <= 1; j++ )
			{
      				unsigned char k = input[row+i][col+j];
      				accumulation += k * filter[1+i][1+j];
      				weightsum += filter[1+i][1+j];
    			}
  		}
  		output[row][col] = (unsigned char)(accumulation/weightsum);
	}
  }
  return output;
}


int main(int argc, char **argv)
{
	SDoublePlane image1 =  SImageIO::read_png_file("cat.png");
	SDoublePlane image2 =  SImageIO::read_png_file("tiger.png");

	cout<<"Read the images\n";
	SDoublePlane real1;
	SDoublePlane real2;
	
	SDoublePlane imaginary1;
	SDoublePlane imaginary2;
	
	//fft(image1,real1,imaginary1);
	//fft(image2,real2,imaginary2);

	

	SDoublePlane gaussian_filter(3,3);
	gaussian_filter[0][0] = 1/16.0;
	gaussian_filter[0][1] = 1/8.0;
	gaussian_filter[0][2] = 1/16.0;
	gaussian_filter[1][0] = 1/8.0;
	gaussian_filter[1][1] = 1/4.0;
	gaussian_filter[1][2] = 1/8.0;
	gaussian_filter[2][0] = 1/16.0;
	gaussian_filter[2][1] = 1/8.0;
	gaussian_filter[2][2] = 1/16.0;
	
	cout<<"Built the filter\n";

	SDoublePlane lowPass = convolve_general(image1, gaussian_filter); //Applying gaussian filter, we get the image after low pass filtered

	SDoublePlane highPass = convolve_general(image2, gaussian_filter);

	cout<<"Applied Filters\n";
	for(int i=0;i<image2.rows();i++)
	{
		for(int j=0;j<image2.cols();j++)
		{
			highPass[i][j] = image2[i][j] - highPass[i][j];
		}
	}

	SDoublePlane merged(lowPass.rows(),lowPass.cols());
	
	for(int i=0;i<merged.rows();i++)
	{
		for(int j=0;j<merged.cols();j++)
		{
			merged[i][j] = (lowPass[i][j]+highPass[i][j])/2;
		}
	}
	SImageIO::write_png_file("merged.png",merged,merged,merged);
}








