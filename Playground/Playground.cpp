// Playground.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
using namespace std;

/**********************************************************************/
/**************************** Helper function *************************/
/**********************************************************************/

vector<double>  MatrixVectorMultiplication(vector<vector<double>> matrix, vector<double> fieldVector)
{
    vector<double> resultVector(matrix.size(), 0.);


    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
        {
            resultVector[j] += fieldVector[j] * matrix[i][j];
        }
        rotate(fieldVector.begin(), fieldVector.begin() + 1, fieldVector.end());
    }
    return resultVector;
}


/*
* Creates a matrix with alternating sign and same value. Except that there are two more positive then negative values per row. Thus the matrix-vektor multiplikation
*/
void CreateMatrix(vector<vector<double>>& matrix)
{
    if (matrix.size() % 2 != 0)
    {
        throw logic_error("Scheme requires a sqaure matrix with even dimension.");
    }
    double val = 0.55;
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
        {
            if (j == i)
            {
                matrix[i][j] = val;
            }
            else
            {
                matrix[i][j] = val * pow(-1, j + i + 1);
            }
        }
    }
}

template <typename T>
inline vector<T> get_diagonal(int position, vector<vector<T>> U)
{

    vector<T> diagonal(U.size());

    int k = 0;
    // U(0,l) , U(1,l+1), ... ,  U(n-l-1, n-1)
    for (int i = 0, j = position; (i < U.size() - position) && (j < U.size()); i++, j++)
    {
        diagonal[k] = U[i][j];
        k++;
    }
    for (int i = U.size() - position, j = 0; (i < U.size()) && (j < position); i++, j++)
    {
        diagonal[k] = U[i][j];
        k++;
    }

    return diagonal;
}


int main()
{
    int slot_count = 100;
    vector<vector<double>> curlDummy(slot_count, vector<double>(slot_count));
    vector<vector<double>> matrixDia(slot_count, vector<double>(slot_count));
    vector<double> fieldDummy(slot_count);

    CreateMatrix(curlDummy);

    for (int i = 0; i < slot_count; i++)
    {
        matrixDia[i] = get_diagonal(i, curlDummy);
    }



    //Pay attention that the final value doesn't grow too much. If scaleVal is set to 1.5 the accuracy drops after 8 itterations, 
    //since the 20 bits of pre dot resolution is not sufficient enough anymore.
    double scaleVal = 1.1;
    double value = 2.01234567890123456789;




    for (int i = 0; i < slot_count; i++)
    {
        fieldDummy[i] = value;
    }

    fieldDummy = MatrixVectorMultiplication(matrixDia, fieldDummy);
}
