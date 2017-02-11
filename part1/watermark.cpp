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
#include <fstream>

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
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag)
{
	
}
// Normalizes into [0-255]
void normalize(SDoublePlane &input)
{
	double min,max;
	int count = 0;
	min = max = input[0][0];
	for(int i=0;i<input.rows();i++)
	{
		for(int j=0;j<input.cols();j++)
		{
			if(input[i][j]<0.0)
				count++;
			if(min>input[i][j])
				min = input[i][j];
			else if(max<input[i][j])
				max = input[i][j];
		}
	}

	for(int i=0;i<input.rows();i++)
	{
		for(int j=0;j<input.cols();j++)
		{
			input[i][j] = ((input[i][j]-min)*255)/(max-min);
		}
	}
}
// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input);


// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);


int main(int argc, char **argv)
{
	try
	{
		if(argc < 4)
		{
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    p2 problemID inputfile outputfile" << endl;
			return -1;
      		}
		string part = argv[1];
    		string inputFile = argv[2];
    		string outputFile = argv[3];
    		cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
	
		SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    		if(part == "1.1")
		{
			SDoublePlane fft_real(input_image.rows(), input_image.cols());
			SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
			SDoublePlane E(input_image.rows(), input_image.cols());
			fft(input_image,fft_real,fft_imaginary);

			cout<< input_image.rows() << ' ' << input_image.cols() << endl;

			for(int row = 0;row<input_image.rows();row++)
			{
				for(int col=0;col<input_image.cols();col++)
				{
					E[row][col] = log10(sqrt(pow(fft_real[row][col],2.0f) + pow(fft_imaginary[row][col],2.0f)));
				}
			}
			normalize(E);
			SImageIO::write_png_file(outputFile.c_str(),E,E,E);
      		}
    		else if(part == "1.2")
      		{
			SDoublePlane fft_real(input_image.rows(), input_image.cols());
			SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
			SDoublePlane E(input_image.rows(), input_image.cols());
			fft(input_image,fft_real,fft_imaginary);

			cout<< input_image.rows() << ' ' << input_image.cols() << endl;

			for(int row=352;row<=356;row++)
			{
				for(int col=0;col<input_image.cols();col++)
				{
					fft_real[row][col] = fft_imaginary[row][col] = 0.0;
				}
			}

			for(int row=156;row<=160;row++)
			{
				for(int col=0;col<input_image.cols();col++)
				{
					fft_real[row][col] = fft_imaginary[row][col] = 0.0;
				}
			}

			SDoublePlane output_image(input_image.rows(), input_image.cols());
			ifft(fft_real,fft_imaginary,output_image);
			SImageIO::write_png_file(outputFile.c_str(),output_image,output_image,output_image);
      		}
    		else if(part == "1.3")
      		{
			if(argc < 6)
	  		{
				cout << "Need 6 parameters for watermark part:" << endl;
	    			cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    			return -1;
	  		}
			SDoublePlane fft_real(input_image.rows(), input_image.cols());
			SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
			fft(input_image,fft_real,fft_imaginary);

			string op(argv[4]);
			string s(argv[5]);
			int N;
			sscanf(s.c_str(),"%d",&N);
			srand(N);
			if(op == "add")
	  		{
	    			// add watermark
				int l = 40;
				double* v = new double[40];
				for(int i=0;i<l;i++)
				{
					v[i] = (rand()%2);
				}
				for(int row=0;row<fft_real.rows();row++)
				{
					for(int col=0;col<fft_real.cols();col++)
					{
						fft_real[row][col] = fft_real[row][col] + alpha*abs(fft_real[row][col])*
					}
				}

				SDoublePlane output_image(input_image.rows(), input_image.cols());
				ifft(fft_real,fft_imaginary,output_image);
				SImageIO::write_png_file(outputFile.c_str(),output_image,output_image,output_image);
	  		}
			else if(op == "check")
	  		{
	    			// check watermark
	  		}
			else
	  			throw string("Bad operation!");
			int N = atoi(argv[5]);
		}
		else
			throw string("Bad part!");
	}
	catch(const string &err)
	{
    		cerr << "Error: " << err << endl;
  	}
}