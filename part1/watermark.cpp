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

double PI = 3.14159;
int l = 20;
double r = 110.0;
double alpha = 10.0;
double t = -0.5;

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

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag)
{
	SDoublePlane E(fft_real.rows(), fft_real.cols());
	for(int row = 0;row<fft_real.rows();row++)
	{
		for(int col=0;col<fft_real.cols();col++)
		{
			E[row][col] = 1+log(sqrt(pow(fft_real[row][col],2) + pow(fft_imag[row][col],2)));
		}
	}
	normalize(E);
	return E;
}
// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input_image)
{
	SDoublePlane fft_real(input_image.rows(), input_image.cols());
	SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
	SDoublePlane E(input_image.rows(), input_image.cols());
	fft(input_image,fft_real,fft_imaginary);

	cout<< input_image.rows() << ' ' << input_image.cols() << endl;

	/*ofstream ofs;
  	ofs.open("fft_real.csv",ofstream::out);
	for(int row=0;row<fft_real.rows();row++)
	{
		ofs<<fft_real[0][0];
		for(int col=0;col<fft_real.cols();col++)
		{
			ofs<<","<<fft_real[row][col];
		}
		ofs<<endl;
	}
	ofs<<endl;
	ofs.close();*/

	SDoublePlane before(input_image.rows(), input_image.cols());
	before = fft_magnitude(fft_real,fft_imaginary);
	normalize(before);
	SImageIO::write_png_file("before.png",before,before,before);

	/*double sum=0.0;
	int count=0;
	for(int row=0;row<fft_real.rows();row++)
	{
		for(int col=0;col<fft_real.cols();col++)
		{
			sum += fft_real[row][col];
			count++;
		}
	}
	double newFreq = sum/count;
	cout<<"New value: "<<newFreq<<endl;*/

	for(int row=352;row<=356;row++)
	{
		for(int col=0;col<input_image.cols();col++)
		{
			fft_real[row][col] = fft_real[351][col];
			fft_imaginary[row][col] = fft_imaginary[351][col];
			//fft_real[row][col] = newFreq ;
			//fft_imaginary[row][col] = 0.0;
		}
	}

	for(int row=156;row<=160;row++)
	{
		for(int col=0;col<input_image.cols();col++)
		{
			fft_real[row][col] = fft_real[155][col];
			fft_imaginary[row][col] = fft_imaginary[155][col];
			//fft_real[row][col] = newFreq;
			//fft_imaginary[row][col] = 0.0;
		}
	}

	SDoublePlane output_image(input_image.rows(), input_image.cols());
	ifft(fft_real,fft_imaginary,output_image);

	SDoublePlane after(input_image.rows(), input_image.cols());
	after = fft_magnitude(fft_real,fft_imaginary);
	normalize(after);
	SImageIO::write_png_file("after.png",after,after,after);

	return output_image;
}

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input_image, int N)
{
	srand(N);
	double* v = new double[l];
	for(int i=0;i<l;i++)
	{
		v[i] = (rand()%2);
	}
	int cx = input_image.rows()/2;
	int cy = input_image.cols()/2;
	int d;

	SDoublePlane fft_real(input_image.rows(), input_image.cols());
	SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
	fft(input_image,fft_real,fft_imaginary);
	
	//double part = r/(l/2);
	//int x,y;
	for(int i=0;i<l/2;i++)
    	{
		//y = part*i;
		//x = (int)sqrt(r*r-y*y);
        	float theta=(float)(i)*(360.0/l);
		int x= (int)(r* cos((theta*PI)/180));
        	int y= (int)(r* sin((theta*PI)/180));
		fft_real[cx+x][cy+y] = fft_real[cx+x][cy+y] + alpha*fabs(fft_real[cx+x][cy+y])*v[i];
		fft_real[cx-x][cy-y] = fft_real[cx-x][cy-y] + alpha*fabs(fft_real[cx-x][cy-y])*v[i];
		
		fft_real[cx+x][cy-y] = fft_real[cx+x][cy-y] + alpha*fabs(fft_real[cx+x][cy-y])*v[i+(l/2)];
		fft_real[cx-x][cy+y] = fft_real[cx-x][cy+y] + alpha*fabs(fft_real[cx-x][cy+y])*v[i+(l/2)];

		//fft_real[cx+x][cy+y] = fft_real[cx-x][cy-y] = fft_real[cx+x][cy-y] = fft_real[cx-x][cy+y] = 20.0;
	}
	
	/*for(int k=0;k<(l/2);k++)
	{
		d = (int)floor(sqrt(r*r - 4*k*k));
		fft_real[cx-d][cy+2*k] = fft_real[cx-d][cy+2*k] + alpha*abs(fft_real[cx-d][cy+2*k])*v[k];
		fft_real[cx+d][cy-2*k] = fft_real[cx+d][cy-2*k] + alpha*abs(fft_real[cx+d][cy-2*k])*v[k];
	}

	for(int k=0;k<(l/2);k++)
	{
		d = (int)floor(sqrt(r*r-4*k*k));
		fft_real[cx-d][cy-2*k] = fft_real[cx-d][cy-2*k] + alpha*abs(fft_real[cx-d][cy-2*k])*v[k+(l/2)];
		fft_real[cx+d][cy+2*k] = fft_real[cx+d][cy+2*k] + alpha*abs(fft_real[cx+d][cy+2*k])*v[k+(l/2)];
	}*/

	SDoublePlane after(input_image.rows(), input_image.cols());
	after = fft_magnitude(fft_real,fft_imaginary);
	normalize(after);
	SImageIO::write_png_file("spectWatermarked.png",after,after,after);
	
	SDoublePlane output_image(input_image.rows(), input_image.cols());
	ifft(fft_real,fft_imaginary,output_image);
	return output_image;
}

double getMean(double array[], int size)
{
         double sum = 0.0;
         for(int i = 0; i < size; i++)
         {
             sum += array[i];
         }

        return (sum / size);
}

double getStDDev(double array[], int size)
{
    double mean = getMean(array,size);
    
    double sum = 0.0;
    for(int i = 0; i < size; i++)
    {
        sum += (array[i]-mean)*(array[i]-mean);
    }
    return sqrt(sum / size);
}

double getCovariance(double a[],double b[],int size)
{
    double mean1 = getMean(a,size);
    double mean2 = getMean(b,size);
    double sum=0;
    
    for(int i = 0; i < size; i++)
    {
        sum += (a[i]-mean1)*(b[i]-mean2);
    }
    return (sum/size);
}


double getCorrelation(double a[], double b[],int size)
{
    return (getCovariance(a,b,size)/(getStDDev(a,size)*getStDDev(b,size)));
}

// Write this in Part 1.3 -- check if watermark N is in image
bool check_image(const SDoublePlane &input_image, int N)
{
	SDoublePlane fft_real(input_image.rows(), input_image.cols());
	SDoublePlane fft_imaginary(input_image.rows(), input_image.cols());
	fft(input_image,fft_real,fft_imaginary);
	
	double* v = new double[l];
    	double* c = new double[l];
	srand(N);
	for(int i=0;i<l;i++)
    	{
        	v[i] = (rand()%2);
	}
	int cx = input_image.rows()/2;
    	int cy = input_image.cols()/2;
    
    	for(int i=0;i<l/2;i++)
	{
        	float theta = (float)(i)*(360.0/l);
       		int x = r* cos((theta*PI)/180);
        	int y = r* sin((theta*PI)/180);

        	c[i] = fft_real[cx+x][cy+y];
        	c[i+(l/2)] = fft_real[cx-x][cy+y];
	}
	double correlation = getCorrelation(c,v,l);
	cout<<correlation<<endl;
	if(correlation>=t)
    		return true;
    	else
    		return false;
}


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

			E = fft_magnitude(fft_real,fft_imaginary);

			double sum = 0.0;
			int count=0;
			for(int i=0;i<E.rows();i++)
			{
				for(int j=0;j<E.cols();j++)
				{
					sum += E[i][j];
					count++;
				}
			}
			cout<<sum/count<<endl;
			SImageIO::write_png_file(outputFile.c_str(),E,E,E);
      		}
    		else if(part == "1.2")
      		{
			SDoublePlane output_image(input_image.rows(), input_image.cols());
			output_image = remove_interference(input_image);
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
			/*SDoublePlane sample(input_image.rows(), input_image.cols());
			for(int row=0;row<input_image.rows();row++)
			{
				for(int col=0;col<input_image.cols();col++)
				{
					sample[row][col] = 10;
				}
			}
			int cx = sample.rows()/2;
			int cy = sample.cols()/2;
			double PI = 3.14159;
			int l = 10;
			double r = 110;
			double alpha = 10.0;
			double t=0.4;

			sample[cx][cy] = 240;
			for(int i=0;i<l/2;i++)
    			{
        			float theta=(float)(i+1)*(360.0/l);
				int x= (int)(r* cos((theta*PI)/180));
        			int y= (int)(r* sin((theta*PI)/180));
				//cout<<"Setting: "<<(cx+x)<<" "<<cy+y<<endl;
				sample[cx+x][cy+y] = 240;//fft_real[cx+x][cy+y] + alpha*abs(fft_real[cx+x][cy+y])*v[i];
				//cout<<"Setting: "<<(cx-x)<<" "<<(cy-y)<<endl;
				sample[cx-x][cy-y] = 240;//fft_real[cx-x][cy-y] + alpha*abs(fft_real[cx-x][cy-y])*v[i];
		
				//cout<<"Setting: "<<(cx+x)<<" "<<cy-y<<endl;
				sample[cx+x][cy-y] = 240;//fft_real[cx+x][cy-y] + alpha*abs(fft_real[cx+x][cy-y])*v[i+(l/2)];
				//cout<<"Setting: "<<(cx-x)<<" "<<cy+y<<endl;
				sample[cx-x][cy+y] = 240;//fft_real[cx-x][cy+y] + alpha*abs(fft_real[cx-x][cy+y])*v[i+(l/2)];
			}
			
			SImageIO::write_png_file("sample.png",sample,sample,sample);*/
			string op(argv[4]);
			int N = atoi(argv[5]);
			if(op == "add")
	  		{
	    			// add watermark
				SDoublePlane output_image(input_image.rows(), input_image.cols());
				output_image = mark_image(input_image,N);
				SImageIO::write_png_file(outputFile.c_str(),output_image,output_image,output_image);

				/*int fp;
				int total=0;
				srand(time(NULL));
				for(int i=0;i<10;i++) //10 runs
				{
					fp=0;
            				int random = rand();
					output_image = mark_image(input_image,random);
					for(int j=0;j<100;j++)
					{
						if(check_image(output_image,rand()%1000))
                                        		fp++;
					}
            				cout<<"i:"<<i+1<< ", N:"<<random<<", FP:"<<(double)fp/100<<endl;
					total += fp;
				}
				cout<<total/(100.0*10)<<endl;*/
	  		}
			else if(op == "check")
	  		{
	    			// check watermark
				SDoublePlane output_image(input_image.rows(), input_image.cols());
				if(check_image(input_image,N))
				{
					cout<<"Found Watermark"<<endl;
				}
				else
				{
					cout<<"Could not find any watermark"<<endl;
				}
	  		}
			else
	  			throw string("Bad operation!");
		}
		else
			throw string("Bad part!");
	}
	catch(const string &err)
	{
    		cerr << "Error: " << err << endl;
  	}
}