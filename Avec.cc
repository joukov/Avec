/*********************************************************************
Copyright 2000 - 2006 Bruno Wittmer
          2006 Michael Henke
This file is part of the Avec package.

Avec is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Avec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Avec; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**********************************************************************/

#include "Avec.h"
#include "Avec2D.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <cstdlib>


// Initialize the Palette for MultiGraph
// The first entry is the makrker type, the second one the color for markers and lines
std::vector< std::pair < int, int> > Init_palette()
{
	std::vector<int> the_colors;
	the_colors.push_back(2);
	the_colors.push_back(3);
	the_colors.push_back(4);
	the_colors.push_back(6);
	the_colors.push_back(7);
	the_colors.push_back(8);
	the_colors.push_back(9);

	std::vector< std::pair < int, int> > retval;
	for(int i = 0; i < 30; i++){
		retval.push_back(std::pair<int, int>(20+i%9, the_colors[i%the_colors.size()]));
	}
	return retval;
}

int Avec::VERBOSE_FLAG=1;
int Avec::OUTPUT_PRECISION = 8;
int Avec::OUTPUT_WIDTH = 10;
std::vector< std::pair < int, int> > Avec::MGRAPH_PALETTE = Init_palette();


//! Root Class Implementation
/*!
Detailed description
Usage with ROOT:\n
Use SetName("name"); to name your Avec in ROOT.\n
Use Write() to write the Avec to File.\n
You can also use Write("another_name")\n
 */
ClassImp(Avec)
;

////////////////////////////
// Implementation of Avec //
//           ||           //
//           vv           //

//! Range Constructor  
/*!
This Constructor creates an Avec object of size 'newsize' and fills it with values from 'low' to 'high'.\n
Example:
\verbatim
  Avec vec1(4,0.2,1.7);   // Vector with 4 elements filled with 0.2 - 1.7
  vec1.print();
  >> 0.2
  >> 0.7
  >> 1.2
  >> 1.7
\endverbatim
*/
Avec::Avec(size_t newsize, dattype low, dattype high) : std::vector<dattype>(newsize)
{
  dattype step = (high-low)/(newsize-1);
  dattype val = low-step;
  for(iterator i = begin(); i!=end(); i++){
    *i = val += step;
  }
}

//! Constructor for Avec2D conversion
/*!
This Constructor creates an Avec by concatenating all rows of the Avec2D object
*/
Avec::Avec(const Avec2D & source)
{
  for(Avec2D::const_iterator it = source.begin();it != source.end(); it++)
    *this=(*this)&(*it);
  //cout << "Created Avec with " << size() << " entries" << std::endl;
}


//! Constructor with array initialization
Avec::Avec(size_t newsize, const dattype data[]): std::vector<dattype>(newsize)
{
  for(iterator i=begin();i!=end();i++,data++){
    *i = *data;
  }
}

//! Constructor from STL vector of floats
Avec::Avec(const std::vector<float> & source)
{
  std::back_insert_iterator<Avec> ii(*this);

  std::copy(source.begin(), source.end(), ii);

  //for(std::vector<float>::const_iterator i=source.begin();i!=source.end();i++)
  //  push_back(*i);
}


// Useful methods

void Avec::help(const std::string s) // Overview of methods (incomplete!)
{
  if(equal(s.begin(),s.end(), std::string("vecread").begin()))
    std::cout << "vecread(std::vector<dattype> &vec, const char * filename)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("vecwrite").begin()))
    std::cout << "vecwrite(std::vector<dattype> &vec, const char * filename)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("operator+").begin()))
    std::cout << "Avec operator+(const Avec& source)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("operator-").begin()))
    std::cout << "Avec operator-(const Avec& source)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("operator&").begin()))
    std::cout << "Avec operator&(const Avec& source)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("print").begin()))
    std::cout << "void print(ostream &out)" << std::endl;
  if(equal(s.begin(),s.end(), std::string("help").begin()))
    std::cout << "void help(string s)" << std::endl;
}

//! Print contents to a stream
void Avec::print(std::ostream & out, const char* sep) const
{
  copy(begin(),end(), std::ostream_iterator<dattype>(out,sep));
}

std::ostream& operator<< (std::ostream& o, const Avec& v)
{
  v.print(o,"\t");
  return o;
}

 
///////////////
// Operators //
///////////////
//!  Returns the range from first to last.
/*!
 * \param first first index of the range.
 * \param last last index of the range.
 * \return A vector with a copy of the range first:last.\n

 * Example: \n
\verbatim
  Avec vec1(10,1,10);   // Vector with 10 elements filled with 1 - 10
  Avec vec2=vec1(4,7);  // Copy the range 4:7 to vec2
  vec2.print();
  >> 5
  >> 6
  >> 7
  >> 8
\endverbatim
*/
Avec Avec::operator()(size_t first, size_t last) const
{
  if(empty()) return Avec();
  if(first <= 0) first = 0;
  if(first >= size()) first = size() - 1;
  last++;
  if(last <= 0) last = size();
  if(last > size()) last = size();
  if(last < first) return Avec();
  return Avec(begin() + first, begin() + last);
}

/*Avec Avec::operator[](AvecRange r) const
{
  if(empty())return Avec();
  if(r.first <= 0)r.first=0;
  if(r.first >= size())r.first=size()-1;
  r.last++;
  if(r.last<=0)r.last=0;
  if(r.last>size())r.last=size();
  if(r.last<r.first)return Avec();
  Avec result=*this;
  //return Avec(begin()+first,begin()+last);
  return Avec(result.begin()+r.first,result.begin()+r.last);
}*/

//!  Returns the range from first to last.
// void Avec::operator=(const std::vector<float>& source)
// {
//   this->clear();
//   //back_insert_iterator<Avec> ii(*this);
//   insert(this->end(), source.begin(), source.end());
//   return;
// }

// Arithmetic with vectors
Avec Avec::operator+(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f+*g;
  return result;
}

Avec Avec::operator-(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f-*g;
  return result;
}

Avec Avec::operator*(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=(*f)*(*g);
  return result;
}

Avec Avec::operator/(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=(*f)/(*g);
  return result;
}

Avec Avec::operator%(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=(long int)(*f)%(long int)(*g);
  return result;
}

void Avec::operator+=(const Avec& source)
{
  size_t length= std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i+=*g;
}

void Avec::operator-=(const Avec& source)
{
  size_t length= std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i-=*g;
}

void Avec::operator*=(const Avec& source)
{
  size_t length= std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i*=*g;
}

void Avec::operator/=(const Avec& source)
{
  size_t length= std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i/=*g;
}

void Avec::operator%=(const Avec& source)
{
  size_t length= std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i = (long int)*i  % (long int)*g;
}

// Arithmetic with scalars

Avec Avec::operator+(double offset) const
{
  Avec result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i+offset;
  return result;
}

Avec Avec::operator-(double offset) const
{
  Avec result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i-offset;
  return result;
}

Avec Avec::operator*(dattype factor) const
{
  Avec result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i*factor;
  //transform(begin(), end(),result.begin(), bind2nd(multiplies<dattype>(), factor));
  return result;
}

Avec Avec::operator/(double factor) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i/factor;
  return result;
}

Avec Avec::operator%(long int factor) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=(long int)*i%factor;
  return result;
}

Avec operator+(double offset, const Avec& source)
{
  return source+offset;
}

Avec operator-(double offset, const Avec& source)
{
  return source*-1+offset;
}

Avec operator*(double factor, const Avec& source)
{
  return source*factor;
}

Avec operator/(double factor, const Avec& source)
{
  return Avec(source.size(),factor)/source;
}

Avec Avec::operator-() const
{
  return *this*-1;
}

void Avec::operator+=(double offset)
{
  for(iterator i=begin();i!=end();i++)
    *i+=offset;
}

void Avec::operator-=(double offset)
{
  for(iterator i=begin();i!=end();i++)
    *i-=offset;
}

void Avec::operator*=(double factor)
{
  for(iterator i=begin();i!=end();i++)
    *i*=factor;
}

void Avec::operator/=(double factor)
{
  for(iterator i=begin();i!=end();i++)
    *i/=factor;
}

void Avec::operator%=(long int factor)
{
  for(iterator i=begin();i!=end();i++)
    *i = (long int)*i % factor;
}

// Comparisons



//! Return an Avec of same size filled with 1 and 0 depending on the comparison for each element
Avec Avec::operator<(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f<*g;
  return result;
}

Avec Avec::operator>(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f>*g;
  return result;
}

Avec Avec::operator<=(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f<=*g;
  return result;
}

Avec Avec::operator>=(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f>=*g;
  return result;
}

Avec Avec::operator<(dattype value) const
{
  Avec result(size(),0.0);
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i<value) *f=1;
  //transform(begin(), end(),result.begin(), bind2nd(less<dattype>(), value));
  return result;
}

Avec Avec::operator>(dattype value) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i>value) *f=1;
  return result;
}

Avec Avec::operator<=(dattype value) const
{
  Avec result(size(),0.0);
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i<=value) *f=1;
  //transform(begin(), end(),result.begin(), bind2nd(less<dattype>(), value));
  return result;
}

Avec Avec::operator>=(dattype value) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i>=value) *f=1;
  return result;
}

Avec Avec::operator!=(dattype value) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i!=value) *f=1;
  return result;
}

Avec Avec::operator==(dattype value) const
{
  Avec result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    if(*i==value) *f=1;
  return result;
}

Avec Avec::operator&&(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f&&*g;
  return result;
}

Avec Avec::operator||(const Avec& source) const
{
  size_t length= std::min(size(),source.size());
  Avec result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i = *f || *g;
  return result;
}



// Other operators


//           AA           //
//           ||           //
// Implementation of Avec //
////////////////////////////

//////////////////////
// Helping functions//
//////////////////////


//! Concatenate vectors
Avec operator&(const Avec& vec1, const Avec& vec2)
{
  Avec result=vec1;
  copy(vec2.begin() , vec2.end() , back_inserter( result ) ) ;
  return result;
}


//! Create list of vectors. e.g. vecwrite((v1,v2,v3,v4),"outfile.dat") 

std::vector<Avec*> operator|(Avec& first,Avec& second)
{
  std::vector<Avec*> result;
  result.push_back(&first);
  result.push_back(&second);
  return result;
}

std::vector<Avec*> operator|(Avec& first, const std::vector<Avec*>& second)
{
  std::vector<Avec*> result=second;
  //result.push_back(&first);
  result.insert(result.begin(),&first);
  return result;
}

std::vector<Avec*> operator|(const std::vector<Avec*>& first, Avec& second)
{
  std::vector<Avec*> result=first;
  result.push_back(&second);
  return result;
}

//! Remove entries that are nan
void avec_remove_nan(Avec& data, dattype replacement)
{
  for(Avec::iterator it = data.begin(); it != data.end(); it++) if(*it != *it) *it = replacement;
}

//! Return an Avec that contains all entries of data, for which mask is non-zero
Avec avec_mask(const Avec& data, const Avec& mask)
{
  Avec retval;

  if(data.size() != mask.size()){
    std::cerr << "Error in avec_mask: data and mask do not have the same size, data.size(): " << data.size() << " ; mask.size(): " << mask.size() << std::endl;
    return retval;
  }
  for(Avec::size_type i = 0; i < data.size(); i++) if(mask[i] != 0) retval.push_back(data[i]);
  return retval;
}



//! Write vector to a file
void vecwrite(const std::vector<dattype> &vec, const std::string &filename)
{
  if(Avec::VERBOSE_FLAG) std::cout << "Writing vector to file " << filename << std::endl;
  fstream file(filename.c_str(), std::ios::out);
  if(!file)return;
  file << std::setprecision(Avec::OUTPUT_PRECISION);
  copy(vec.begin(),vec.end(), std::ostream_iterator<dattype>(file,"\n"));
  file.close();
}

//! Write vector list to a stream
std::ostream& operator<< (std::ostream& o, const std::vector<Avec*>& list)
{
  if(!list.empty()){

    // Determine smallest vector in list
    size_t min=list[0]->size();
    for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++)
      if((*it)->size()<min)min=(*it)->size();
    
    // Write list into file
    o << std::setprecision(Avec::OUTPUT_PRECISION);
    for(size_t f=0; f<min;f++){
      for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++)
	//o << setprecision(8) << setw(10) << (*(*it))[f] << "  ";
	o << (*(*it))[f] << "  ";
      o << std::endl;
    }
  }
  return o;
}

//! Write vector list to a file
void vecwrite(const std::vector<Avec*> &list, const std::string &filename)
{
  if(list.empty())return;

  if(Avec::VERBOSE_FLAG) std::cout << "Writing vector list to file " << filename << std::endl;
  fstream file(filename.c_str(), std::ios::out);
  if(!file)return;

  file << list;
  file.close();
}


//! Read vector from a file
int vecread(std::vector<dattype> &vec, const std::string &filename)
{
  fstream file(filename.c_str(), std::ios::in);
  if(!file){
    std::cerr << "Could not open file " << filename << std::endl;
    return -1;
  }
  vec.clear();
  dattype buffer=0;
  while(file >> buffer){
    vec.push_back(buffer);
  }
  if(file.eof()){
    if(Avec::VERBOSE_FLAG) 
      std::cout << "Vector created with size " << vec.size() << std::endl;
  }
  else
    std::cerr << "Error while reading file " << filename << std::endl;
  file.close();
  return 0;
}

//! Read vector list from a file
int vecread(const std::vector<Avec*> &list, const std::string &filename)
{
  if(list.empty()){
    std::cerr << "Empty list!" << std::endl;
    return -1;
  }
  fstream file(filename.c_str(), std::ios::in);
  if(!file){
    std::cerr << "Could not open file " << filename << std::endl;
    return -1;
  }

  // Determine smallest vector in list
  size_t min=list[0]->size();
  for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++)
    if((*it)->size()<min)min=(*it)->size();
  //cout << "Reading " << min << " Elements" << std::endl;

  if(min != 0){
    for(size_t f=0; f<min;f++){
      for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++){
	file >> (*(*it))[f];
	//cout << (*(*it))[f] << "\n";
      }
    }
  }
  else{
    for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++)
      (*it)->clear();
    double buff = 0;
    int ready = 0;
    while (!ready){
      for(std::vector<Avec*>::const_iterator it=list.begin();it!=list.end();it++){
	if(file >> buff){
	  (*it)->push_back(buff);
	}
	else{
	  ready = 1;
	  break;
	}
      }
    }
    if(Avec::VERBOSE_FLAG) std::cout << "Created vectors of size " << (*list.begin())->size() << std::endl;
  }
  file.close();
  return 0;
}


//! Sum of all entries
dattype vsum(const std::vector<dattype> &vec)
{
  //dattype result=0;
  //for(Avec::const_iterator i=vec.begin(); i!=vec.end(); i++)
  //  result+=*i;
  //return result;
  return accumulate(vec.begin(),vec.end(),(dattype)0);
}

//! Maximum Element
dattype vmax(const std::vector<dattype> &vec)
{
  if(vec.empty())return 0;
  return *max_element(vec.begin(),vec.end());
}

//! Minimum Element
dattype vmin(const std::vector<dattype> &vec)
{
  if(vec.empty())return 0;
  return *min_element(vec.begin(),vec.end());
}

//! Index (Location) of maximum element
size_t lvmax(const std::vector<dattype> &vec)
{
  return max_element(vec.begin(),vec.end())-vec.begin();
}

//! Index (Location) of minimum element
size_t lvmin(const std::vector<dattype> &vec)
{
  return std::min_element(vec.begin(),vec.end())-vec.begin();
}


//! function to "differentiate" an Avec. 
/*! In the output Avec one finds in each position (except the first) the 
   difference between this position and the previous one of the input Avec
*/
Avec vdiff(const Avec& v1)
{
  Avec v2 = v1;
  for( Avec::size_type i = 1;i < v1.size(); i++ )
    {
      v2[i] = v1[i] - v1[i-1];
    }
  if(!v2.empty()) v2[0] = 0;
  return v2;
}

Avec vdiff2(const Avec& v1)
{
  Avec::size_type sz=v1.size();
  Avec v2(sz);
  if( sz > 1 ){
    v2[0]=v1[1]-v1[0];
    v2[sz-1]=v1[sz-1]-v1[sz-2];
    for( Avec::size_type i = 1;i < sz-1; i++ )
      {
	v2[i] = v1[i+1] - v1[i-1];
      }
  }
  
  return v2;
}

//! function to "integrate" an Avec. 
/*! In the output Avec one finds in each position the sum of 
   all positions from the first up to this position of the input Avec.
*/
Avec vint(const Avec& v1)
{
  Avec v2 = v1;
  for( Avec::size_type i = 1;i < v1.size(); i++ )
    {
      v2[i] = v2[i-1] + v1[i];
    }
  return v2;
}

//! function to calculate the mean of an Avec
dattype vmean(const Avec& v1)
{
  dattype mean = 0.0;
  for( Avec::size_type i = 0; i < v1.size(); i++){
      mean += v1[i]*i;
  }
  mean /=  vsum(v1);
  return mean;
}

//! function to calculate the RMS of an Avec
dattype vrms(const Avec& v1)
{
  return sqrt(vsum(v1*v1)/v1.size());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Compress the Avec pair to a smaller size with regular x spacing
/*!
A regular grid of n x-intervals is created and all entries of the original Avec are assigned to one of them
All original entries assigned to one new entry are averaged (x and y).
If there are gaps in the original x-values that create gaps in the new range, the new empty entries can be removed by setting the flag remove_empty to true
 */
void avec_compress(Avec& x, Avec& y, Avec::size_type n, bool remove_empty)
{
  if(x.size() != y.size()){
    std::cerr << "Error in avec_compress: x.size = " << x.size() << " and y.size = " << y.size() << "  are not the same" << std::endl;
    return;
  }

  if( n >= x.size()){
    std::cerr << "Error in avec_compress: n = " << n << " larger or equal than x.size = " << x.size() << std::endl;
    return;
  }

  double xmin = vmin(x);
  double xmax = vmax(x);
  double dx = (xmax - xmin) / n;

  Avec xsum(n);
  Avec ysum(n);
  Avec norm(n);

  for(Avec::size_type i = 0; i < x.size(); i++){
    Avec::size_type idx = std::min( (Avec::size_type)((x[i] - xmin) / dx), n-1);
    xsum[idx] += x[i];
    ysum[idx] += y[i];
    norm[idx]++;
  }

  xsum /= norm;
  ysum /= norm;

  if(remove_empty){
    x = avec_mask(xsum, norm != 0);
    y = avec_mask(ysum, norm != 0);
  }
  else{
    x = xsum;
    y = ysum;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Compress the Avec pair and the errors to a smaller size with regular x spacing
/*!
A regular grid of n x-intervals is created and all entries of the original Avec are assigned to one of them
All original entries assigned to one new entry are averaged (x and y).
The errors are defined as the rms of the averaging process.
If there are gaps in the original x-values that create gaps in the new range, the new empty entries can be removed by setting the flag remove_empty to true
 */
void avec_compress(Avec& x, Avec& y, Avec& err_y, Avec::size_type n, bool remove_empty)
{
  if(x.size() != y.size()){
    std::cerr << "Error in avec_compress: x.size = " << x.size() << " and y.size = " << y.size() << "  are not the same" << std::endl;
    return;
  }

  if(x.size() != err_y.size()){
    std::cerr << "Error in avec_compress: x.size = " << x.size() << " and err_y.size = " << err_y.size() << "  are not the same" << std::endl;
    return;
  }

  if( n >= x.size()){
    std::cerr << "Error in avec_compress: n = " << n << " larger or equal than x.size = " << x.size() << std::endl;
    return;
  }

  double xmin = vmin(x);
  double xmax = vmax(x);
  double dx = (xmax - xmin) / n;

  Avec xsum(n);
  Avec ysum(n);
  Avec y2sum(n);
  Avec norm(n);

  for(Avec::size_type i = 0; i < x.size(); i++){
    Avec::size_type idx = std::min( (Avec::size_type)((x[i] - xmin) / dx), n-1);
    xsum[idx] += x[i];
    ysum[idx] += y[i];
    y2sum[idx] += y[i]*y[i];
    norm[idx]++;
  }

  y2sum = sqrt((norm*y2sum - ysum*ysum)/norm/norm/(norm-1));
  xsum /= norm;
  ysum /= norm;

  if(remove_empty){
    x = avec_mask(xsum, norm != 0);
    y = avec_mask(ysum, norm != 0);
    err_y = avec_mask(y2sum, norm != 0);
  }
  else{
    x = xsum;
    y = ysum;
    err_y = y2sum;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Compress the Avecs to a smaller size with regular x spacing
/*!
A regular grid of n x-intervals is created and all entries of the original Avecs are assigned to one of them
All original entries assigned to one new entry are averaged (x and y).
If there are gaps in the original x-values that create gaps in the new range, the new empty entries can be removed by setting the flag remove_empty to true
 */
void avec_compress(Avec& x, std::vector<Avec*> yv, Avec::size_type n, bool remove_empty)
{
  std::cout << "Reducing Avecs of size " << x.size() << " to " << n << std::endl;

  Avec xwork;
  for(unsigned int i = 0; i < yv.size(); i++){
    xwork = x;
    avec_compress(xwork, *(yv[i]), n, remove_empty);
  }
  x = xwork;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Compress the Avecs and the errors to a smaller size with regular x spacing
/*!
A regular grid of n x-intervals is created and all entries of the original Avec are assigned to one of them
All original entries assigned to one new entry are averaged (x and y).
The errors are defined as the rms of the averaging process.
If there are gaps in the original x-values that create gaps in the new range, the new empty entries can be removed by setting the flag remove_empty to true
 */
void avec_compress(Avec& x, std::vector<Avec*> yv, std::vector<Avec*> err_yv, Avec::size_type n, bool remove_empty)
{
  if(yv.size() != err_yv.size()){
    std::cerr << "Error in void avec_compress(Avec&, std::vector<Avec*>, std::vector<Avec*>, Avec::size_type, bool): "
	      << "yv.size()=" << yv.size()
	      << " does not match err_yv.size()=" << err_yv.size() << std::endl;
    return;
  }
  std::cout << "Reducing Avecs of size " << x.size() << " to " << n << std::endl;

  Avec xwork;
  for(unsigned int i = 0; i < yv.size(); i++){
    xwork = x;
    avec_compress(xwork, *(yv[i]), *(err_yv[i]), n, remove_empty);
  }
  x = xwork;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! Resize the Avec
WINAVEC_API Avec avec_resize(const Avec& v, Avec::size_type nsize)
{
  Avec result;
  if(v.empty() || nsize < 1)return result;
  double factor=(double)nsize/(double)v.size();
  //cout << "Factor: " << factor << std::endl;

  result=Avec(nsize);
  //cout << "New Size: " << result.size() << std::endl;

  Avec::size_type iold=0;
  Avec::size_type inew=0;

  while(iold < v.size() && inew < nsize){
    //cout << "iold: " << iold << " inew: " << inew << std::endl;
    dattype diff=(inew + 1)/factor-(iold+1);
    //cout << "diff: " << diff << std::endl;
    if(diff >=0){
      result[inew]+=v[iold];
      iold++;
    }
    else{
      result[inew]+=v[iold]*(diff+1);
      inew++;
      if(inew /factor>iold && inew < nsize){
	result[inew]+=v[iold]*-diff;
	iold++;
      }
    }
  }

  result *= factor;
  return result;

}

//! Resize the Avec
WINAVEC_API Avec rescale_size(const Avec& v, double factor)
{
  if(v.empty())return Avec();
  Avec::size_type n=v.size();
  if(Avec::VERBOSE_FLAG) std::cout << "Vector size: " << n << std::endl;
  n=(int)((double)n*factor);
  if(Avec::VERBOSE_FLAG) std::cout << "New Vector size: " << n << std::endl;

  return avec_resize(v,n);
}

//! Create an Avec with y-values, x1 should be sorted!!!
// x1 and y1 contain a function, x2 contains a set of x-values, which should be covered by x1
// The returned Avec contains y-values corresponding to the x2-values, interpolating linearly between pairs (x1,y1)
Avec avec_map(const Avec& x1, const Avec& y1, const Avec& x2)
{
  if(x1.empty() || x1.size() != y1.size()) return Avec(x2.size());

  Avec y2(x2.size());

  for(Avec::size_type i = 0; i < x2.size(); i++){
    if( x2[i] < x1[0]){
      y2[i] = y1[0];
      continue;
    }
    if( x2[i] >= x1.back()){
      y2[i] = y1.back();
      continue;
    }

    for(Avec::size_type f = 0; f < x1.size() - 1 ; f++){
      if(x2[i] >= x1[f] && x2[i] < x1[f+1]){
	y2[i] = y1[f] + (y1[f+1] - y1[f])/(x1[f+1] - x1[f])*(x2[i] - x1[f]);
      }
    }
  }
  return y2;
}




//! math.h functions
/*!
\verbatim
Avec acos(const std::vector<dattype> &vec)
Avec asin(const std::vector<dattype> &vec)
Avec atan(const std::vector<dattype> &vec)
Avec atan2(const std::vector<dattype> &y, const std::vector<dattype> &x)
Avec ceil(const std::vector<dattype> &vec)
Avec cos(const std::vector<dattype> &vec)
Avec cosh(const std::vector<dattype> &vec)
Avec exp(const std::vector<dattype> &vec)
Avec fabs(const std::vector<dattype> &vec)
Avec floor(const std::vector<dattype> &vec)
Avec fmod(const std::vector<dattype> &x, dattype y)
Avec log(const std::vector<dattype> &vec)
Avec log10(const std::vector<dattype> &vec)
Avec pow(const std::vector<dattype> &vec, dattype y)
Avec sin(const std::vector<dattype> &vec)
Avec sinh(const std::vector<dattype> &vec)
Avec sqrt(const std::vector<dattype> &vec)
Avec tan(const std::vector<dattype> &vec)
Avec tanh(const std::vector<dattype> &vec)
\endverbatim
*/
Avec fabs(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=fabs(*i);
  return result;
}

Avec acos(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=acos(*i);
  return result;
}

Avec asin(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=asin(*i);
  return result;
}

Avec atan(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=atan(*i);
  return result;
}

Avec atan2(const std::vector<dattype> &y, const std::vector<dattype> &x)
{
  Avec::size_type rsize=0;
  if(y.size() != x.size())
    std::cerr << "atan2 function got vectors with different size!" << std::endl;
  else
    rsize=y.size();

  Avec result(rsize);

  Avec::iterator f=result.begin();
  Avec::const_iterator g=x.begin();
  for(Avec::const_iterator i=y.begin();i!=y.end();i++,f++,g++)
    *f=atan2(*i,*g);
  return result;
}

Avec ceil(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=ceil(*i);
  return result;
}

Avec cos(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=cos(*i);
  return result;
}

Avec cosh(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=cosh(*i);
  return result;
}

Avec exp(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=exp(*i);
  return result;
}

Avec floor(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=floor(*i);
  return result;
}

Avec fmod(const std::vector<dattype> &x, dattype y)
{
  Avec result(x.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=x.begin();i!=x.end();i++,f++)
    *f=fmod(*i,y);
  return result;
}

Avec log(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=log(*i);
  return result;
}

Avec log10(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=log10(*i);
  return result;
}

Avec pow(const std::vector<dattype> &vec, dattype y)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=pow(*i,y);
  return result;
}

Avec sin(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=sin(*i);
  return result;
}

Avec sinh(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=sinh(*i);
  return result;
}

Avec sqrt(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=sqrt(*i);
  return result;
}

Avec tan(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=tan(*i);
  return result;
}

Avec tanh(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=tanh(*i);
  return result;
}

//////////////////////////
//    Random Numbers    //
//////////////////////////

//! Random numbers from 0 to RAND_MAX
WINAVEC_API Avec rand(std::vector<dattype>::size_type size)
{
  Avec result(size);

  Avec::iterator f=result.begin();
  for(f=result.begin();f!=result.end();f++)
    *f=rand();

  return result;
}

//! Random numbers from 0 to RAND_MAX (obsolete, use rand(std::vector<dattype>::size_type size))
Avec rand(const std::vector<dattype> &vec)
{
  return rand(vec.size());
}

//! Random numbers from 0 to 1
Avec urand(const std::vector<dattype> &vec)
{
  Avec result(vec.size());

  Avec::iterator f=result.begin();
  for(Avec::const_iterator i=vec.begin();i!=vec.end();i++,f++)
    *f=(double)rand()/(double)RAND_MAX;
  return result;
}

//! Random Numbers with gaussian distribution (Box-Muller)
Avec grand(const std::vector<dattype> &vec, double sigma)
{
	//std::cout << "Called Avec grand(const std::vector<dattype> &vec, double sigma)" << std::endl;
	return grand(vec.size(), sigma);
}

//! Random Numbers with gaussian distribution (Box-Muller)
Avec grand(Avec::size_type size, double sigma)
{
	//std::cout << "Called Avec grand(Avec::size_type size, double sigma)" << std::endl;
  Avec result(size);
  for(Avec::iterator f=result.begin();f!=result.end();f++)
    *f=  cos((double)rand()/(double)RAND_MAX*2*M_PI)*sqrt(log((double)rand()/(double)RAND_MAX)*-2)*sigma;
  return result;
}


/////////////////////////////////////////////////////////////////////////

#ifdef __AVECROOT__
// Root-specific stuff
#include "TArrayD.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TLegend.h"
#include "TFile.h"

// Overloaded TObject methods
void Avec::Draw(Option_t* option)
{
  gPad->Clear();
  avec_draw(*this);
  gPad->Update();
}

void Avec::Browse(TBrowser* b)
{
  std::string canvas_name="avec_draw_cv";
  // Check if TCanvas object with this name already exists, if yes remove it
  TCanvas* g_test = (TCanvas*)(gDirectory->FindObject( canvas_name.c_str() ));
  if( g_test != 0 ){
    if(Avec::VERBOSE_FLAG) std::cout << "Reusing TCanvas" << std::endl;
    g_test->cd();
  }
  else {
    if(Avec::VERBOSE_FLAG) std::cout << "creating new TCanvas" << std::endl;
    new TCanvas(canvas_name.c_str(),canvas_name.c_str());
  }
  // Draw the Avec
  Draw(b ? b->GetDrawOption() : "");
}


Avec empty_avec;

//! Retrieve an Avec from a root file through TFile reference
Avec& avec_get(const std::string& name, TFile& file)
{
  Avec* Avec_ptr = 0;
  file.GetObject(name.c_str(),Avec_ptr);
  if(!Avec_ptr){
    std::cerr << "Did not find " << name << std::endl;
    return empty_avec;
  }
  return *Avec_ptr;
}

//! Retrieve an Avec from a root file through file name
Avec& avec_get(const std::string& vecname, const std::string& filename)
{
  TFile f(filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return empty_avec;
  }
  return avec_get(vecname, f);
}



Avec TAD2AV(TArrayD & source)
{
  Avec result(source.GetSize());

  for (unsigned int f=0;f<result.size();f++)
    result[f]=source.At(f);
  return result;
}

//##########################################################################################
//# root function to build a vector of the Avec class out of a TH1 histogram               #
//##########################################################################################

Avec TH12avec(TH1& h1)
{
  Avec vdata(0);    // (0) needed for valid initialisation
  
  int nstrips = h1.GetNbinsX();
     
  double bin_content = 0.0;
  for(int i=0 ; i<nstrips ; i++){
    bin_content = h1.GetBinContent(i+1);   // histo bin index starts at 1!
    vdata.push_back(bin_content);          // vector index starts at 0!
  } 

  return vdata;
}
      
//#################################################################################################
//# root function to create a TH1D histo out of an avec vector                                    #
//#################################################################################################

TH1D* avec2TH1D(Avec& v1, const std::string& histo_title, const std::string histo_name = "avec2histo1")
{
      if(Avec::VERBOSE_FLAG) std::cout << "histo_name.c_str() is:    " << histo_name.c_str() << std::endl;
      
      int nbins = v1.size();

      TH1D* h_test = (TH1D*)(gDirectory->FindObject( histo_name.c_str() ));
      if( h_test != 0 )
      {
          h_test->Delete();
          if(Avec::VERBOSE_FLAG) std::cout << "existing histo deleted!" << std::endl;
      }

      TH1D* h1 = new TH1D( histo_name.c_str() , histo_title.c_str() , nbins , 0.5 , nbins + 0.5 );
 
      for(int i=0; i<nbins; i++)
      {
          h1->SetBinContent( i+1 , v1[i]);
      }
 
      return h1; 
}

//-------------------------------------------------------------------------------------------------

TGraph* create_graph(const Avec& x, const Avec& y)
{
    Avec::size_type size_x = x.size();
    Avec::size_type size_y = y.size();

    if( size_x != size_y){
      std::cerr << "Error in create_graph: Avecs have not the same size" << std::endl;
      std::cerr << "size_x: " << size_x << "  size_y: " << size_y << std::endl;
      return 0;
    }
    if(size_x == 0)return 0; //Skip empty entries
    
	TGraph *gr = new TGraph(size_x);
    // Fill in the data
    for(Avec::size_type i=0; i != size_x; i++) gr->SetPoint(i, x[i], y[i]);
    //gr->Sort();
    return gr;
}

TGraphErrors* create_error_graph(const Avec& x, const Avec& y, const Avec& ex, const Avec& ey)
{
    Avec::size_type size_x = x.size();
    Avec::size_type size_y = y.size();
    Avec::size_type size_ex = ex.size();
    Avec::size_type size_ey = ey.size();

    if( size_x != size_y || size_x != size_ey || size_x != size_ex){
      std::cerr << "Error in create_graph: Avecs have not the same size" << std::endl;
      std::cerr << "size_x: " << size_x << "  size_y: " << size_y << "  size_ex: " << size_ex << "  size_ey: " << size_ey << std::endl;
      return 0;
    }
    if(size_x == 0)return 0; //Skip empty entries
    // Create new TGraphErrors object;
    TGraphErrors *gr = new TGraphErrors(size_x);
    // Fill in the data
    for(Avec::size_type i=0; i < size_x; i++){
    	gr->SetPoint(i, x[i], y[i]);
    	gr->SetPointError(i, ex[i], ey[i]);
    }
    //gr->Sort();
    return gr;
}

//*****************************************************************
//* Helping function to get axis labels and drawing options right *
//*****************************************************************
void avec_draw_paint(TMultiGraph* mgr, const std::string& xtitle, const std::string& ytitle, const std::string& options)
{
  if(!mgr)return;
  if( mgr->GetListOfGraphs()->IsEmpty() ) return;

  // The following lines are needed, or else ROOT crashes
  // For some reason TMultiGraph needs to be drawn before setting axis titles 
  // This is because the axes are only created when the multigraph is drawn.

  // Another Problem: std::string.c_str() creates only a temporary buffer, ROOT doesn't seem to handle that properly for the axis titles and options
  // So we make a persistent copy. (and create a memory leak...)

  std::string & xtitle_copy = * new std::string(xtitle);
  std::string & ytitle_copy = * new std::string(ytitle);
  std::string & options_copy = * new std::string(options);

  mgr->Draw(options_copy.c_str());
  mgr->GetXaxis()->SetTitle(xtitle_copy.c_str());
  mgr->GetYaxis()->SetTitle(ytitle_copy.c_str());
  gPad->Modified();
  return;
}


//#############################################################################
//# root function to display a single Avec vector in a TGraph object          #
//#############################################################################

TGraph *avec_draw(const Avec& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string y_tit = yTitle == "__NO_AXIS_TITLE__" ? y.GetName() : yTitle;
  return avec_draw(Avec(y.size(),1,y.size()), y, title, xTitle, y_tit, options, draw, name);
}

//#############################################################################
//# root function to display two Avec vectors in a TGraph object              #
//#############################################################################

TGraph *avec_draw(const Avec& x, const Avec& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  // Check if vectors are of same size      
  Avec::size_type n = x.size();
  if (n != y.size()){
    std::cerr << "The vectors have to be of same size!" << std::endl;
    return 0;
  }

  std::string graph_name = name != "" ? name : (std::string)y.GetName() + "%" + x.GetName();
  graph_name = GetValidObjectName<TGraph>(graph_name);
  if(Avec::VERBOSE_FLAG) std::cout << "TGraph name is:    " << graph_name << std::endl;

  std::string y_tit = yTitle == AVEC::no_axis_title ? y.GetName() : yTitle;
  std::string x_tit = xTitle == AVEC::no_axis_title ? x.GetName() : xTitle;

  // Create new TGraph object;
  TGraph *gr = new TGraph(n);
  gr->SetName(graph_name.c_str());
  gr->SetTitle(title.c_str());
  // Fill in the data
  for(Avec::size_type i=0; i<n; i++)
    gr->SetPoint(i, x[i], y[i]);
  //gr->Sort();
  gr->SetMarkerStyle(21);
  gr->SetDrawOption(options.c_str());
  gr->GetXaxis()->SetTitle(x_tit.c_str());
  gr->GetYaxis()->SetTitle(y_tit.c_str());

  if(draw) gr->Draw(options.c_str());

  return gr;
}

//###################################################################################
//# root function to display two Avec vectors in a TGraphErrors object              #
//###################################################################################

TGraphErrors *avec_draw(const Avec& x, const Avec& y, const Avec& e_x, const Avec& e_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string graph_name = name != "" ? name : (std::string)y.GetName() + "%" + x.GetName();
  graph_name = GetValidObjectName<TGraph>(graph_name);
  if(Avec::VERBOSE_FLAG) std::cout << "TGraphErrors name is:    " << graph_name << std::endl;

  // Check if vectors are of same size      
  Avec::size_type n = x.size();
  if (n != y.size() || n != e_x.size() || n != e_y.size() ){
    std::cerr << "The vectors have to be of same size!" << std::endl;
    return 0;
  }

  // Create new TGraphErrors object;
  TGraphErrors *gr = new TGraphErrors(n);
  gr->SetName(graph_name.c_str());
  gr->SetTitle(title.c_str());
  // Fill in the data
  for(Avec::size_type i=0; i<n; i++){
    gr->SetPoint(i, x[i], y[i]);
    gr->SetPointError(i, e_x[i], e_y[i]);
  }
  //gr->Sort();
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetTitle(xTitle.c_str());
  gr->GetYaxis()->SetTitle(yTitle.c_str());

  if(draw) gr->Draw(options.c_str());

  return gr;
}


//###########################################################################################
//# root function to display a set of Avec vectors against one Avec in a TMultiGraph object #
//###########################################################################################
TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec*>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::vector<const Avec*> ny;
  for(std::vector<Avec*>::size_type i=0; i<y.size(); i++)
    ny.push_back(y[i]);
  return avec_draw(x, ny, title, xTitle, yTitle, options, draw, name);
}

//###########################################################################################
//# root function to display a set of Avec vectors against one Avec in a TMultiGraph object #
//###########################################################################################
TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::vector<const Avec*> ny;
  for(std::vector<Avec>::size_type i=0; i<y.size(); i++)
    ny.push_back(&y[i]);
  return avec_draw(x, ny, title, xTitle, yTitle, options, draw, name);
}


//###########################################################################################
//# root function to display a set of Avec vectors against one Avec in a TMultiGraph object #
//###########################################################################################
TMultiGraph *avec_draw(const Avec& x, const std::vector<const Avec*>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string graph_name = GetValidObjectName<TMultiGraph>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TMultiGraph name is:    " << graph_name.c_str() << std::endl;

  // Check if vectors are of same size      
  //Avec::size_type n = x.size();

  // Create new TMultiGraph object;
  TMultiGraph *mgr = new TMultiGraph();
  mgr->SetName(graph_name.c_str());
  mgr->SetTitle(title.c_str());

  for(std::vector<Avec*>::size_type it=0;it!=y.size();it++){
    Avec::size_type n = y[it]->size();
    if (n != x.size()){
      std::cout << "The vectors have to be of same size!" << std::endl;
      return 0;
    }

    TGraph *gr = new TGraph(n);
    //gr->SetName(name.c_str());
    //gr->SetTitle(title.c_str());
    // Fill in the data
    for(Avec::size_type i=0; i != n; i++)
      gr->SetPoint(i, x[i], (*y[it])[i]);
    //gr->Sort();
    //gr->SetMarkerStyle(21+(it%5));
    //gr->SetMarkerColor(1+(it+1)%9);
    //gr->SetLineColor(1+(it+1)%9);
    if(Avec::MGRAPH_PALETTE.empty()){
    	gr->SetMarkerStyle(21+(it%5));
    	gr->SetMarkerColor(1+(it+1)%9);
    	gr->SetLineColor(1+(it+1)%9);
    }
    else{
    	std::pair<int, int>& pal_entry = Avec::MGRAPH_PALETTE[it%Avec::MGRAPH_PALETTE.size()];
    	gr->SetMarkerStyle(pal_entry.first);
    	gr->SetMarkerColor(pal_entry.second);
    	gr->SetLineColor(pal_entry.second);
    }
    mgr->Add(gr);
  }

  if(draw)avec_draw_paint(mgr, xTitle, yTitle, options);

  return mgr;
}


//#######################################################################################################
//# root function to display a set of Avec vectors with errors against one Avec in a TMultiGraph object #
//#######################################################################################################
TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec*>& y, const std::vector<Avec*>& e_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  if(y.size() != e_y.size()){
    std::cerr << "Error in avec_draw: The sets of y-Avecs and of error Avecs have not the same size" << std::endl;
    return 0;
  }

  std::vector<const Avec*> ny, ne_y;
  for(std::vector<Avec*>::size_type i=0; i<y.size(); i++){
    ny.push_back(y[i]);
    ne_y.push_back(e_y[i]);
  }
  return avec_draw(x, ny, ne_y, title, xTitle, yTitle, options, draw, name);
}

//#######################################################################################################
//# root function to display a set of Avec vectors with errors against one Avec in a TMultiGraph object #
//#######################################################################################################
TMultiGraph *avec_draw(const Avec& x, const std::vector<const Avec*>& y, const std::vector<const Avec*>& e_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string graph_name = GetValidObjectName<TMultiGraph>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TMultiGraph name is:    " << graph_name.c_str() << std::endl;

  // Check if vector sets are of same size
  if (y.size() != e_y.size()){
    std::cerr << "Error in avec_draw: The sets of y-Avecs and of error Avecs have not the same size" << std::endl;
    return 0;
  }

  // Create new TMultiGraph object;
  TMultiGraph *mgr = new TMultiGraph(graph_name.c_str(),title.c_str());
  //mgr->SetName(name.c_str());
  //mgr->SetTitle(title.c_str());

  for(std::vector<Avec*>::size_type it=0;it!=y.size();it++){
    Avec::size_type n = y[it]->size();
    TGraphErrors *gr = new TGraphErrors(n);
    //gr->SetName(name.c_str());
    //gr->SetTitle(title.c_str());
    // Fill in the data
    for(Avec::size_type i=0; i != n; i++){
      gr->SetPoint(i, x[i], (*y[it])[i]);
      gr->SetPointError(i, 0, (*e_y[it])[i]);
    }
    //gr->Sort();
    //gr->SetMarkerStyle(21+(it%5));
    //gr->SetMarkerColor(1+(it+1)%9);
    //gr->SetLineColor(1+(it+1)%9);
    if(Avec::MGRAPH_PALETTE.empty()){
    	gr->SetMarkerStyle(21+(it%5));
    	gr->SetMarkerColor(1+(it+1)%9);
    	gr->SetLineColor(1+(it+1)%9);
    }
    else{
    	std::pair<int, int>& pal_entry = Avec::MGRAPH_PALETTE[it%Avec::MGRAPH_PALETTE.size()];
    	gr->SetMarkerStyle(pal_entry.first);
    	gr->SetMarkerColor(pal_entry.second);
    	gr->SetLineColor(pal_entry.second);
    }
    
    mgr->Add(gr);
  }

  if(draw)avec_draw_paint(mgr, xTitle, yTitle, options);

  return mgr;
}

//#######################################################################################################
//# root function to display a set of Avec vectors against another set of Avecs in a TMultiGraph object #
//#######################################################################################################
TMultiGraph *avec_draw(const std::vector<Avec*>& x, const std::vector<Avec*>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  // Check if sets are of same size      
  if (x.size() != y.size()){
    std::cerr << "Error in avec_draw: The sets x and y have not the same size!" << std::endl;
    return 0;
  }

  std::vector<const Avec*> nx, ny;
  for(std::vector<Avec*>::size_type i=0; i<x.size(); i++){
    nx.push_back(x[i]);
    ny.push_back(y[i]);
  }
  return avec_draw(nx, ny, title, xTitle, yTitle, options, draw, name);

}

//#######################################################################################################
//# root function to display a set of Avec vectors against another set of Avecs in a TMultiGraph object #
//#######################################################################################################
TMultiGraph *avec_draw(const std::vector<const Avec*>& x, const std::vector<const Avec*>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string graph_name = GetValidObjectName<TMultiGraph>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TMultiGraph name is:    " << graph_name.c_str() << std::endl;

  // Check if sets are of same size      
  if (x.size() != y.size()){
    std::cerr << "Error in avec_draw: The sets x and y have not the same size!" << std::endl;
    return 0;
  }

  // Create new TMultiGraph object;
  TMultiGraph *mgr = new TMultiGraph();
  mgr->SetName(graph_name.c_str());
  mgr->SetTitle(title.c_str());

  for(std::vector<Avec*>::size_type idx=0; idx!=y.size(); idx++){
    Avec::size_type size_x = x[idx]->size();
    Avec::size_type size_y = y[idx]->size();
    if( size_x != size_y){
      std::cerr << "Error in avec_draw: Avecs Nr." << idx << "have not the same size and will be skipped" << std::endl;
      std::cerr << "size_x: " << size_x << "  size_y: " << size_y << std::endl;
      continue;
    }
    if( size_x == 0){
      if(Avec::VERBOSE_FLAG > 0) std::cerr << "Warning in avec_draw: Avecs Nr." << idx << "are empty and will be skipped" << std::endl;
      continue;
    }
    TGraph *gr = new TGraph(size_x);
    //gr->SetName(name.c_str());
    //gr->SetTitle(title.c_str());
    // Fill in the data
    for(Avec::size_type i=0; i != size_x; i++)
      gr->SetPoint(i, (*x[idx])[i], (*y[idx])[i]);
    //gr->Sort();
    //gr->SetMarkerStyle(21+(idx%5));
    //gr->SetMarkerColor(1+(idx+1)%9);
    //gr->SetLineColor(1+(idx+1)%9);
    if(Avec::MGRAPH_PALETTE.empty()){
    	gr->SetMarkerStyle(21+(idx%5));
    	gr->SetMarkerColor(1+(idx+1)%9);
    	gr->SetLineColor(1+(idx+1)%9);
    }
    else{
    	std::pair<int, int>& pal_entry = Avec::MGRAPH_PALETTE[idx%Avec::MGRAPH_PALETTE.size()];
    	gr->SetMarkerStyle(pal_entry.first);
    	gr->SetMarkerColor(pal_entry.second);
    	gr->SetLineColor(pal_entry.second);
    }
    mgr->Add(gr);
  }

  if(draw)avec_draw_paint(mgr, xTitle, yTitle, options);

  return mgr;
}

//#######################################################################################################
//# root function to display a set of Avec vectors against another set of Avecs in a TMultiGraph object #
//#######################################################################################################
TMultiGraph *avec_draw(const std::vector<Avec>& x, const std::vector<Avec>& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  std::string graph_name = GetValidObjectName<TMultiGraph>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TMultiGraph name is:    " << graph_name.c_str() << std::endl;

  // Check if sets are of same size      
  if (x.size() != y.size()){
    std::cerr << "Error in avec_draw: The sets x and y have not the same size!" << std::endl;
    return 0;
  }

  // Create new TMultiGraph object;
  TMultiGraph *mgr = new TMultiGraph();
  mgr->SetName(graph_name.c_str());
  mgr->SetTitle(title.c_str());

  for(std::vector<Avec>::size_type idx=0; idx!=y.size(); idx++){
    Avec::size_type size_x = x[idx].size();
    Avec::size_type size_y = y[idx].size();
    if( size_x != size_y){
      std::cerr << "Error in avec_draw: Avecs Nr." << idx << " have not the same size and will be skipped" << std::endl;
      std::cerr << "size_x: " << size_x << "  size_y: " << size_y << std::endl;
      continue;
    }
    if( size_x == 0){ //Skip empty entries
      if(Avec::VERBOSE_FLAG > 0) std::cerr << "Warning in avec_draw: Avecs Nr." << idx << " are empty and will be skipped" << std::endl;
      continue;
    }
    TGraph *gr = new TGraph(size_x);
    //gr->SetName(name.c_str());
    //gr->SetTitle(title.c_str());
    // Fill in the data
    for(Avec::size_type i=0; i != size_x; i++)
      gr->SetPoint(i, (x[idx])[i], (y[idx])[i]);
    //gr->Sort();
    //gr->SetMarkerStyle(21+(idx%5));
    //gr->SetMarkerColor(1+(idx+1)%9);
    //gr->SetLineColor(1+(idx+1)%9);
    if(Avec::MGRAPH_PALETTE.empty()){
    	gr->SetMarkerStyle(21+(idx%5));
    	gr->SetMarkerColor(1+(idx+1)%9);
    	gr->SetLineColor(1+(idx+1)%9);
    }
    else{
    	std::pair<int, int>& pal_entry = Avec::MGRAPH_PALETTE[idx%Avec::MGRAPH_PALETTE.size()];
    	gr->SetMarkerStyle(pal_entry.first);
    	gr->SetMarkerColor(pal_entry.second);
    	gr->SetLineColor(pal_entry.second);
    }
    mgr->Add(gr);
  }
  
  if(! mgr->GetListOfGraphs() ){
    if(Avec::VERBOSE_FLAG > 0) std::cerr << "Warning in avec_draw: Multigraph is empty and will not be drawn" << std::endl;
  	return mgr;
  }
  
  if( draw && mgr->GetListOfGraphs() )avec_draw_paint(mgr, xTitle, yTitle, options);

  return mgr;
}

//###################################################################################################################
//# root function to display a set of Avec vectors with errors against another set of Avecs in a TMultiGraph object #
//###################################################################################################################
TMultiGraph *avec_draw(const std::vector<Avec>& x, const std::vector<Avec>& y, const std::vector<Avec>& ey, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  // Check if sets are of same size      
  if (x.size() != y.size() || x.size() != ey.size()){
    std::cerr << "Error in avec_draw: The sets x,y and ey have not the same size!" << std::endl;
    return 0;
  }

  std::string graph_name = GetValidObjectName<TMultiGraph>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TMultiGraph name is:    " << graph_name.c_str() << std::endl;

  // Create new TMultiGraph object;
  TMultiGraph *mgr = new TMultiGraph();
  mgr->SetName(graph_name.c_str());
  mgr->SetTitle(title.c_str());

  for(std::vector<Avec>::size_type idx=0; idx!=y.size(); idx++){
  	TGraphErrors* gr = create_error_graph(x[idx], y[idx], Avec(x[idx].size()), ey[idx]);
  	if(!gr) continue;
    //gr->SetMarkerStyle(21+(idx%5));
    //gr->SetMarkerColor(1+(idx+1)%9);
    //gr->SetLineColor(1+(idx+1)%9);
     if(Avec::MGRAPH_PALETTE.empty()){
    	gr->SetMarkerStyle(21+(idx%5));
    	gr->SetMarkerColor(1+(idx+1)%9);
    	gr->SetLineColor(1+(idx+1)%9);
    }
    else{
    	std::pair<int, int>& pal_entry = Avec::MGRAPH_PALETTE[idx%Avec::MGRAPH_PALETTE.size()];
    	gr->SetMarkerStyle(pal_entry.first);
    	gr->SetMarkerColor(pal_entry.second);
    	gr->SetLineColor(pal_entry.second);
    }
    mgr->Add(gr);
  }

  if(draw)avec_draw_paint(mgr, xTitle, yTitle, options);

  return mgr;
}

// Fill the conetents of an Avec into a histogram and plot it
TH1* avec_plot(const Avec& v, int nbins, const std::string& title, const std::string& name)
{
	std::string hist_name = GetValidObjectName<TH1D>(name);
	if(Avec::VERBOSE_FLAG) std::cout << "TH1D name is:    " << hist_name.c_str() << std::endl;

	// Create the limits for the x-axis
	double xmax=vmax(v);
	double xmin=vmin(v);
	double diff=xmax-xmin;
	xmax+= 0.1*diff;
	xmin-= 0.1*diff;
	if(Avec::VERBOSE_FLAG > 2) std::cout << "Histogram range: " << xmin << " - " << xmax << "  nbins: " << nbins << std::endl;
	
	TH1 * h=new TH1D(hist_name.c_str(),title.c_str(),nbins,xmin,xmax);

	hfill(v,*h);

	h->Draw();
	return h;
}

// Fill the conetents of an Avec into a histogram and plot it, fill only if mask nonzero
TH1* avec_plot_mask(const Avec& v, const Avec& mask, int nbins, const std::string& title, const std::string& name)
{
  if(v.size() != mask.size()){
    std::cerr << "avec_plot_mask: The two Avecs should have the same size!\n" 
	      << "v.size(): " << v.size() << "   mask.size(): " << mask.size() << std::endl;
    return 0;
  }

  std::string hist_name = GetValidObjectName<TH1D>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TH1D name is:    " << hist_name.c_str() << std::endl;

  Avec v2;
  for(Avec::size_type i=0; i<v.size(); i++){
    if(mask[i])v2.push_back(v[i]);
  }

  // Create the limits for the x-axis
  double xmax=vmax(v2);
  double xmin=vmin(v2);
  double diff=xmax-xmin;
  xmax+= 0.1*diff;
  xmin-= 0.1*diff;


  TH1 * h=new TH1D(hist_name.c_str(),title.c_str(),nbins,xmin,xmax);
  hfill(v2, *h);

  h->Draw();
  return h;
}

//! Fill the values of one Avec into a histogram weighting with a second Avec
/*! If w is an empty Avec default weight is 1
*/
void hfill(const Avec& v, TH1& h, const Avec& w)
{
  dattype weight=1;
  int w_flag = (w.size()!=0);
  if(v.size() != w.size() &&   w_flag){
    std::cerr << "Error in hfill(const Avec&, TH1&, const Avec&): The two Avecs should have the same size!\n";
    return;
  }
  for(Avec::size_type i=0; i<v.size(); i++){
    if(w_flag)weight=w[i];
    h.Fill(v[i],weight);
  }
}

//! Fill the values of one Avec into a histogram weighting with a second Avec (OBSOLETE)
void hfill(const Avec& v, const Avec& w, TH1& h)
{
  if(v.size() != w.size()){
    std::cerr << "The two Avecs should have the same size!\n";
    return;
  }
  for(Avec::size_type i=0; i<v.size(); i++)
    h.Fill(v[i],w[i]);
}


// Fill the values of two Avecs into a 2D histogram if mask is not empty ignore values where mask == 0
void hfill(const Avec& x, const Avec& y, TH2& h, const Avec& mask)
{
  //if(x.size() != y.size()){
  //  std::cerr << "Error in hfill(const Avec& x, const Avec& y, TH2& h, const Avec& mask): x.size(): " << x.size() << " y.size(): " << y.size() << " should be the same\n";
  //  return;
  //}
  bool mask_flag = true;
  if(mask.empty())
  	mask_flag = false;
  else{
  	if(mask.size() < std::min(x.size(), y.size())){
	    std::cerr << "Error in hfill(const Avec& x, const Avec& y, TH2& h, const Avec& mask): Avec mask is not of same size as  x and y\n";
	    return;
  	}
  }
  
  Avec::const_iterator i,f, m;
  for(i=x.begin(),f=y.begin(), m=mask.begin(); i!=x.end() && f!=y.end() ;i++,f++,m++)
    if(mask_flag && !*m)continue;
    else h.Fill(*i,*f);
}

// Fill the values of two Avecs into a 2D histogram and draw it
TH2* avec_plot(const Avec& x, const Avec& y, int nbinsx, int nbinsy, const std::string& title, const std::string& name)
{

  // Create the limits for the x-axis
  double xmax=vmax(x);
  double xmin=vmin(x);
  double diffx=xmax-xmin;
  xmax+= 0.1*diffx;
  xmin-= 0.1*diffx;

  // Create the limits for the y-axis
  double ymax=vmax(y);
  double ymin=vmin(y);
  double diffy=ymax-ymin;
  ymax+= 0.1*diffy;
  ymin-= 0.1*diffy;

  std::string hist_name = GetValidObjectName<TH2D>(name);
  TH2 * h=new TH2D(hist_name.c_str(), title.c_str(), nbinsx, xmin, xmax, nbinsy, ymin, ymax);

  //for(Avec::const_iterator i=v.begin();i!=v.end();i++)
  //  h->Fill(*i);
  hfill(x,y,*h);

  h->Draw();
  return h;
}

  // Create the limits for an axis
void avec_get_axis_range(const Avec& data, double& low, double& up)
{
  up = vmax(data);
  low = vmin(data);
  double diff = up - low;
  up+= 0.1*diff;
  low-= 0.1*diff;
  return;
}

WINAVEC_API TH2* avec_plot_mask(const Avec& x, const Avec& y, const Avec& mask, int nbinsx, int nbinsy, const std::string & title, const std::string & name)
{
  // Create the limits for the x-axis
  double xmin, xmax;
  avec_get_axis_range(x, xmin, xmax);

  // Create the limits for the y-axis
  double ymin, ymax;
  avec_get_axis_range(y, ymin, ymax);

  std::string hist_name = GetValidObjectName<TH2D>(name);
  TH2 * h=new TH2D(hist_name.c_str(), title.c_str(), nbinsx, xmin, xmax, nbinsy, ymin, ymax);

  hfill(x, y, *h, mask);

  h->Draw();
  return h;
}

#endif
