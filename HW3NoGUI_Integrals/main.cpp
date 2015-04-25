#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
using namespace std;

/**
 * Will use the Polynomial class to represent a polynomial for the purposes of this problem.
 * 
 * 
 */



/**
 * @class Polynomial
 * @author Nikola
 * @date 04/13/15
 * @file main.cpp
 * @brief Class representing an arbitrary polynomial with real coefficients.
 */
class Polynomial
{
public:
    /**
     * @brief Default constructor. Constant 0 polynomial.
     * @return Polynomial object instance.
     */
    Polynomial()
    {
        data.assign(1,0);
        indeces.assign(1,0);
    }
	/** Constructor using an array of double values.
	    */
	Polynomial(double* a, int num)
    {
        if(a!=NULL && num>0)
        {
            data.assign(num, 0);
            indeces.assign(num, 0);
            for(int i=0;i<num;i++)
            {
                // Filling indeces with positions.
                indeces[i]=i;
                // Copying over the data.
                data[i]=a[i];
            }
        }
    }
	/** Copy constructor. Proper deep copy.
        */
	Polynomial(const Polynomial &other)
    {
        data = other.data;
        indeces = other.indeces;
    }
    Polynomial& operator=(Polynomial other)
    {
        this->data.swap(other.data);
        this->indeces.swap(other.indeces);
        return *this;
    }
    /**
     * @param index - which coefficient to get
     * @return value of index coefficient. n means a_n*x^n.
     */
    double getCoefficient(int index) const
    {
        for(int i=0;i<indeces.size();++i)
            if(indeces[i]==index)
                return data[i];

        // If item not found return 0.
        return 0.0;
    }
    /**
     * @brief Setter.
     * @param index - where to insert. n means a_n*x^n.
     * @param val - Value to insert.
     */
    void setCoefficient(int index, double val)
    {
        //Look through indeces. If found, change. If not, append.
        for(int i=0;i<indeces.size();++i)
            if(indeces[i]==index)
                {data[i]=val; return;}

        indeces.push_back(index);
        data.push_back(val);
    }
    /**
     * @brief Pretty print the polynomial to std::cout.
     */
    void print() const
    {
        for (int i=getDegree();i>=0;--i)
        {
            if(i!=0 && this->getCoefficient(i)!=0)
                cout<<this->getCoefficient(i)<<"*x^"<<i<<" + ";
            if(i==0)
                cout<<this->getCoefficient(i)<<endl;
        }
    }
    Polynomial& operator+=(const Polynomial &other)
    {
        for (int i=0;i<=this->getDegree(other);++i)
            this->setCoefficient(i, this->getCoefficient(i)+other.getCoefficient(i));
        return *this;
    }
    friend Polynomial operator+(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial ret(lhs);
        ret+=rhs;
        return ret;
    }

    Polynomial& operator-=(const Polynomial &other)
    {
        for (int i=0;i<=getDegree(other);++i)
            this->setCoefficient(i, this->getCoefficient(i)-other.getCoefficient(i));
        return *this;
    }
    friend Polynomial operator-(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial ret(lhs);
        ret-=rhs;
        return ret;
    }
    friend Polynomial operator*(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial p;
        for(int i=0;i<=lhs.getDegree(); ++i)
        {
            if(lhs.getCoefficient(i)!=0)
            {
                Polynomial tmp(rhs);
                for(int j=0;j<tmp.indeces.size();++j)
                    tmp.indeces[j]+=i;
                for(int j=0 ; j<tmp.data.size();++j)
                    tmp.data[j]*=lhs.getCoefficient(i);
                p+=tmp;
            }
        }
        return p;
    }
    Polynomial& operator*=(const Polynomial &other)
    {
        (*this)=(*this)*other;
        return *this;
    }
    Polynomial& operator*=(double c)
    {
        for(int i=0;i<data.size();++i)
        {
            data[i]*=c;
        }
        return *this;
    }
    friend Polynomial operator*(double c, const Polynomial& other)
    {
        Polynomial p(other);
        return p*=c;
    }
    friend Polynomial operator*(const Polynomial& other, double c)
    {
        Polynomial p(other);
        return p*=c;
    }
    
    /**
     * @brief Interates a polynomial
     * @return Integral of polynomial with constant set to 0.
     */
    Polynomial& operator++()    // Pre-increment. Integrate
    {
        for(int i=0 ; i<indeces.size() ; i++)
            ++indeces[i];
        for(int i=0;i<data.size();i++)
            data[i]/=indeces[i];
        indeces.push_back(0);
        data.push_back(0);
        return *this;
    }
    
    /**
     * @brief Post-integrates a polynomial.
     * @return Integral of polynomial.
     */
    Polynomial operator++(int) // Post-increment. Integrate
    {
        Polynomial tmp(*this);
        operator++();
        return tmp;
    }
    Polynomial& operator--()   //Pre-decrement. Differentiate
    {
        for(int i=0;i<data.size();++i)
            data[i]*=indeces[i];
        for(int i=0;i<indeces.size();++i)
            --indeces[i];
        for(int i=0;i<indeces.size();++i)
            if(indeces[i]==-1)
            {
                indeces.erase(indeces.begin()+i);
                data.erase(data.begin()+i);
                break;
            }
        return *this;

    }
    Polynomial operator--(int) // Post-decrement. Differentiate
    {
        Polynomial tmp(*this);
        operator--();
        return tmp;
    }
    double operator()(double x) const // Function call operator
    {
        // Implements the Horner scheme for polynomial evaluation.
        int deg = getDegree();
        double p=getCoefficient(deg);
        for(int i=deg-1;i>=0;i--)
        {
            p*=x;
            p+=getCoefficient(i);
        }
        return p;
    }
    double operator[](int index) const //Return index-th coefficient
	{
        return this->getCoefficient(index);
    }

private:
	std::vector<double> data; /**< Vector to hold coefficient data */
    std::vector<int> indeces; /**< Vector to hold index corresp to data coeff */
    int getDegree() const
    {
        int index_max=0;

        // Find largest element in indeces. Corresponds to degree coefficient.
        for(int i=0;i<indeces.size();++i)
            if(indeces[i]>index_max)
                index_max=indeces[i];
        return index_max;
    }
    /**
     * @param other
     * @return Maximum degree of both polynomials.
     */
    int getDegree(const Polynomial &other) const
    {
        return getDegree()> other.getDegree()? getDegree() : other.getDegree();
    }
};


/**
 * End of polynomial class!!!
 * 
 **/

/**
 * @class AlmostPolynomialFunction
 * @author Nikola
 * @date 04/24/15
 * @file main.cpp
 * @brief This will realise the need for a function that is the absolute value integral
 * of a polynomial function.
 */
class AlmostPolynomialFunction
{
public:
    AlmostPolynomialFunction(){Polynomial pol; roots = NULL; n=0;}
    AlmostPolynomialFunction(Polynomial p,double* rts=NULL, int nm=0);
    double operator()(double x) const;
private:
    Polynomial pol;
    double* roots;
    int n;
};

/**
 * @brief Takes a polynomial function and the roots of it sorted in increasing order
 * and the number of roots n.
 * @param p - Polynomial function
 * @param rts - array of double values roots of polynomial sorted in increasing order
 * @param nm - number of roots/length of rts array.
 * @return AlmostPolynomialFunction object
 */
AlmostPolynomialFunction::AlmostPolynomialFunction(Polynomial p,double* rts, int nm)
{
    pol = p; n = nm; 
    if(n==0)
        roots = NULL;
    else
    {
        roots = new double[n]; 
        for(int i=0;i<n;i++) 
            roots[i]=rts[i];
    }
}

/**
 * @brief Evaluates the function at a given x.
 * @param x - value to evaluate at
 * @return Value of the function there.
 */
double AlmostPolynomialFunction::operator()(double x) const
{
    if(roots==NULL)
        return pol(x);
    if(x<=roots[0]+0.000001)
        return pol(x);
    
    double result=0.0; int i=0;
    for(;i<n;i++)
        if(roots[i]<=x+0.000001)
            result+=(i%2==0?2:-2)*pol(roots[i]);
        else
            break;
    result+=(i%2==0?1:-1)*pol(x);
    return result;
}

/**
 * @class Indefinite_Integral
 * @author Nikola
 * @date 04/24/15
 * @file main.cpp
 * @brief Defines an indefinite integral class. Pol is the integrand polynomial.
 */
class Indefinite_Integral
{
public:
    Indefinite_Integral(){pol = *(new Polynomial());}
    Indefinite_Integral(char* str){pol = parseString(str);}
    Polynomial getIntegrand()const{return pol;}
    Polynomial getPrimitive()const{Polynomial ret(pol); return ++ret;} //Sets integration constant to 0.
    AlmostPolynomialFunction getAbsPrimitive(double* sortedroots=NULL, int numroots=0) const;
protected:
    Polynomial pol;
    Polynomial parseString(char* string) const;
};




/**
 * @brief Returns the indefinite integral of a polynomial pol. Will not work properly without
 * a sorted list of the roots of the integrand and their number!
 * @param sortedroots - sorted array of the roots of the integrand. If there are no roots omit argument.
 * e.g. {-1.0, 5.0, 10.0}.
 * @param numroots - number of roots. If none, omit argument.
 * @return AlmostPolynomialFunction object which is the primitive of |P(x)| and evaluates at any x.
 */
AlmostPolynomialFunction Indefinite_Integral::getAbsPrimitive(double* sortedroots, int numroots) const
{
    return AlmostPolynomialFunction(getPrimitive(), sortedroots, numroots);
}

/**
 * @brief Parses a polynomial string representation into a Polynomial object.\
 * expects string in format eg. "5*x^10 +7*x^4 -16*x^2 -20"
 * @param string - char* string of polynomial representation.
 * @return Polynomial object representing the string polynomial.
 */
Polynomial Indefinite_Integral::parseString(char* string) const 
{
    Polynomial ret;
    char* whereami = string;
    while(true)
    {
        double coef = strtod(whereami, &whereami);
        if(*whereami=='\0') {ret.setCoefficient(0,coef); break;}
        whereami+=3;
        int index = (int) strtod(whereami, &whereami);
        ret.setCoefficient(index, coef);
    }
    return ret;
}


class DefiniteIntegral:public Indefinite_Integral
{
public:
    DefiniteIntegral();
    DefiniteIntegral(char* str, double start, double end);
    double getStart()const{return start;}
    double getEnd()const{return end;}
    void setStart(double val){start = val;}
    void setEnd(double val){end = val;}
    double evaluate()const{return getPrimitive()(end) - getPrimitive()(start);}
    double absEvaluate(double* sortedroots=NULL, int numroots=0)const;
protected:
    double start, end;
};

DefiniteIntegral::DefiniteIntegral():Indefinite_Integral()
{
    start = end = 0;
}

DefiniteIntegral::DefiniteIntegral(char* str, double start, double end)\
                                        :Indefinite_Integral(str)
{
    this->start = start; this->end = end;
}

double DefiniteIntegral::absEvaluate(double* sortedroots, int numroots) const
{
    AlmostPolynomialFunction primitive = getAbsPrimitive(sortedroots, numroots);
    return primitive(end) - primitive(start);
}



int main()
{
	

    cout<<endl<<endl<<"**** \t 1. Indefinite integral \t ****"<<endl<<endl;

    

    char initpol[] = "5*x^10 +7*x^4 -16*x^2 -20";
    Indefinite_Integral ii1(initpol);
    ii1.getIntegrand().print();
    ii1.getPrimitive().print();
    DefiniteIntegral di0(initpol, 0, 1);
    cout<<di0.evaluate()<<"(-23.43788)"<<endl<<endl;
    
    char newpol1[] = "1*x^2 -3*x^1 -4"; 
    char newpol2[] = "0.3333*x^3 -1.5*x^2 -4*x^1"; //Integral of the above
    double rts1[] = {-1.0, 4.0};
    DefiniteIntegral di1(newpol1, 0,1);
    AlmostPolynomialFunction apf1(di1.getPrimitive(), rts1, 2);
    di1.getIntegrand().print();
    cout<<"Evaluates the integral of |P(x)| for the function above at -5, 0, 3, 10."<<endl;
    cout<<apf1(-5)<<" "<<apf1(0)<<" "<<apf1(3)<<" "<<apf1(10)<<endl;
    cout<<"Defines and evaluates the definite integral of |P(x)|."<<endl;
    cout<<"Integral from -5 to 3 of 1*x^2 -3*x^1 -4 is:"<<endl;
    DefiniteIntegral di2(newpol1, -5, 3);
    cout<<di2.absEvaluate(rts1, 2)<<"(should be 80)"<<endl;
    
    cout<<endl<<"Indefinite integral of\n\
    (x+5)*(x+4)*(x+3)*(x+2)*(x+1)*x*(x-1)*(x-2)*(x-3)*(x-4)\n\
    which is 1*x^10 +5*x^9 -30*x^8 -150*x^7 +273*x^6 +1365*\n\
x^5 -820*x^4 -4100*x^3 +576*x^2 +2880*x^1"<<endl;
    
    char newpol3[] = "1*x^10 +5*x^9 -30*x^8 -150*x^7 +273*x^6 +1365*x^5 -820*x^4 -4100*x^3 +576*x^2 +2880*x^1";
    DefiniteIntegral di3(newpol3,-4.5,-1.5);
    double roots1[] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4};
    
    cout<<"The integral of this from -4.5 to -1.5 is:"\
    <<di3.absEvaluate(roots1, 10)<<endl;
    cout<<"Should be 13284.8"<<endl;
    di3.getIntegrand().print();
    di3.setStart(-1.5);
    di3.setEnd(2.5);
    cout<<"From -1.5 to 2.5 is: "<<di3.absEvaluate(roots1, 10)<<"(4373.01)"<<endl;
    
    return 0;
}
