#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
using namespace std;

/**
 * Will use the Polynomial class to represent a polynomial for the purposes of this problem.
 * The string will be converted to a polynomial which I will write a function for.
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
	/** Copy constructor.
        */
	Polynomial(const Polynomial &other)
    {
        data = other.data; //Vector in stl copies over data. Deep copy is properly made. Tests below.
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
 * @class Indefinite_Integral
 * @author Nikola
 * @date 04/24/15
 * @file main.cpp
 * @brief Defines an indefinite integral class. Pol is the integrand polynomial.
 */
class Indefinite_Integral
{
public:
    Indefinite_Integral();
    Indefinite_Integral(char* str);
    Polynomial getIntegrand(){return pol;}
    Polynomial getPrimitive(){return pol++;}
    AlmostPolynomialFunction getAbsPrimitive();
protected:
    Polynomial pol;
    Polynomial parseString(char* string) const;
};

Indefinite_Integral::Indefinite_Integral()
{
    pol = *(new Polynomial());
}

Indefinite_Integral::Indefinite_Integral(char* str)
{
    pol = parseString(str);
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
    while(*whereami!='\0')
    {
        double coef = strtod(whereami, &whereami);
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
    DefiniteIntegral(char* str);
    double getStart(){return start;}
    double getEnd(){return end;}
    void setStart(double val){start = val;}
    void setEnd(double val){end = val;}
    double evaluate();
    double absEvaluate();
protected:
    double start, end;
};



int main()
{
	

    cout<<endl<<endl<<"**** \t 1. Indefinite integral ++ -- \t ****"<<endl<<endl;

    

    char initpol[] = "5*x^10 +7*x^4 -16*x^2 -20";
    Indefinite_Integral ara(initpol);
    ara.getIntegrand().print();
    
    
	return 0;
}