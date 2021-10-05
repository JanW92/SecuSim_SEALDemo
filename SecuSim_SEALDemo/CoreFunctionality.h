#pragma once
/*
From the SEALMVDemo project
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <stdlib.h>
#include <vector>
#include "seal/seal.h"
#include "windows.h"
#include "psapi.h"
using namespace std;
using namespace seal;



//Helper function : Prints the parameters in a SEALContext.
inline void print_parameters(const seal::SEALContext & context)
{
    auto& context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

inline string get_parameters(const seal::SEALContext& context)
{
    string parameter;
    auto& context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    parameter= "/\n";
    parameter += "| Encryption parameters :\n" ;
    parameter += "|   scheme: " + scheme_name + "\n";
    parameter += "|   poly_modulus_degree: " + to_string(context_data.parms().poly_modulus_degree()) + "\n" ;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    parameter += "|   coeff_modulus size: ";
    parameter += to_string(context_data.total_coeff_modulus_bit_count()) + " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        parameter += to_string(coeff_modulus[i].bit_count()) + " + ";
    }
    parameter += to_string(coeff_modulus.back().bit_count());
    parameter += ") bits\n";

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        parameter += "|   plain_modulus: " + to_string(context_data.parms().plain_modulus().value()) + "\n";
    }

    parameter += "\\\n";
    return parameter; 
}



// Helper function that prints a matrix (vector of vectors)
template <typename T>
inline void print_full_matrix(vector<vector<T>> matrix, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);
    int row_size = matrix.size();
    int col_size = matrix[0].size();
    for (unsigned int i = 0; i < row_size; i++)
    {
        cout << "[";
        for (unsigned int j = 0; j < col_size - 1; j++)
        {
            cout << matrix[i][j] << ", ";
        }
        cout << matrix[i][col_size - 1];
        cout << "]" << endl;
    }
    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

// Helper function that prints parts of a matrix (only squared matrix)
template <typename T>
inline void print_partial_matrix(vector<vector<T>> matrix, int print_size = 3, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);

    int row_size = matrix.size();
    int col_size = matrix[0].size();

    // Boundary check
    if (row_size < 2 * print_size && col_size < 2 * print_size)
    {
        cerr << "Cannot print matrix with these dimensions: " << to_string(row_size) << "x" << to_string(col_size) << ". Increase the print size" << endl;
        return;
    }
    // print first 4 elements
    for (unsigned int row = 0; row < print_size; row++)
    {
        cout << "\t[";
        for (unsigned int col = 0; col < print_size; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << "..., ";
        for (unsigned int col = col_size - print_size; col < col_size - 1; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << matrix[row][col_size - 1];
        cout << "]" << endl;
    }
    cout << "\t..." << endl;

    for (unsigned int row = row_size - print_size; row < row_size; row++)
    {
        cout << "\t[";
        for (unsigned int col = 0; col < print_size; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << "..., ";
        for (unsigned int col = col_size - print_size; col < col_size - 1; col++)
        {
            cout << matrix[row][col] << ", ";
        }
        cout << matrix[row][col_size - 1];
        cout << "]" << endl;
    }

    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

template <typename T>
inline void print_partial_vector(vector<T> vec, int size, int print_size = 3, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);

    int row_size = size;

    // Boundary check
    if (row_size < 2 * print_size)
    {
        cerr << "Cannot print vector with these dimensions: " << to_string(row_size) << ". Increase the print size" << endl;
        return;
    }

    cout << "\t[";
    for (unsigned int row = 0; row < print_size; row++)
    {
        cout << vec[row] << ", ";
    }
    cout << "..., ";

    for (unsigned int row = row_size - print_size; row < row_size - 1; row++)
    {
        cout << vec[row] << ", ";
    }
    cout << vec[row_size - 1] << "]\n";

    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

template <typename T>
inline void print_full_vector(vector<T> vec, int print_size, int precision = 3)
{
    // save formatting for cout
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(precision);

    cout << "[";
    for (unsigned int x = 0; x < print_size; x++)
    {
        cout << vec[x] << ", ";
    }
    cout << "]" << endl;

    cout << endl;
    // restore old cout formatting
    cout.copyfmt(old_fmt);
}

// Gets a diagonal from a matrix U
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

/*
inline Ciphertext Linear_Transform_Plain(Ciphertext ct, vector<Plaintext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params, bool rescale)
{
    cout << "    + Invoked Linear_Transform_Plain to compute C_vec . P_mat" << endl;

    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;

    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);
    cout << "    + Scale Linear_Transform_Plain add: " << log2(ct_new.scale()) << " bits" << endl;

    vector<Ciphertext> ct_result(U_diagonals.size());
    evaluator.multiply_plain(ct_new, U_diagonals[0], ct_result[0]);

    for (int l = 1; l < U_diagonals.size(); l++)
    {
        Ciphertext temp_rot;
        evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
        evaluator.multiply_plain(temp_rot, U_diagonals[l], ct_result[l]);
    }

    Ciphertext ct_prime;
    evaluator.add_many(ct_result, ct_prime);
    if (rescale)
    {
        cout << "    + Scale Linear_Transform_Plain add_many: " << log2(ct_prime.scale()) << " bits" << endl;
        evaluator.rescale_to_next_inplace(ct_prime);
        cout << "    + after rescale Linear_Transform_Plain add_many: " << log2(ct_prime.scale()) << " bits" << endl;
    }

    return ct_prime;
}
*/
/*
inline Ciphertext Linear_Transform_Cipher(Ciphertext ct, vector<Ciphertext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params, bool rescale, RelinKeys relinKeys)
{
    cout << "    + Invoked Linear_Transform_Cipher C_vec . C_mat:" << endl;
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);

    // Fill ct with duplicate
    Ciphertext ct_rot;

    evaluator.rotate_vector(ct, -U_diagonals.size(), gal_keys, ct_rot);
    Ciphertext ct_new;
    evaluator.add(ct, ct_rot, ct_new);

    cout << "    + Scale add: " << log2(ct_new.scale()) << " bits" << endl;

    vector<Ciphertext> ct_result(U_diagonals.size());
    evaluator.multiply(ct_new, U_diagonals[0], ct_result[0]);

    // Added relin
    if (rescale)
    {
        evaluator.relinearize_inplace(ct_result[0], relinKeys);
    }

    for (int l = 1; l < U_diagonals.size(); l++)
    {
        Ciphertext temp_rot;
        evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
        evaluator.multiply(temp_rot, U_diagonals[l], ct_result[l]);
    }

    Ciphertext ct_prime;
    evaluator.add_many(ct_result, ct_prime);
    if (rescale)
    {
        // Added relin
        evaluator.relinearize_inplace(ct_prime, relinKeys);
        cout << "    + Scale Linear_Transform_Cipher addmany: " << log2(ct_prime.scale()) << " bits" << endl;
        evaluator.rescale_to_next_inplace(ct_prime);
        cout << "    + after rescale Linear_Transform_Cipher addmany: " << log2(ct_prime.scale()) << " bits" << endl;
    }

    return ct_prime;
}
*/
inline void test_Linear_Transformation(int dimension, vector<vector<double>> input_matrix, vector<double> input_vec)
{
    vector<double> result(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            result[k] += input_matrix[i][j] * input_vec[j];
        }
        k++;
    }

    // Print Result vector
    print_partial_vector(result, dimension);
}

//return result
inline vector<double> get_Linear_Transformation_expected_vector(int dimension, vector<vector<double>> input_matrix, vector<double> input_vec)
{
    vector<double> result(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            result[k] += input_matrix[i][j] * input_vec[j];
        }
        k++;
    }

    return result;
}

//return result
inline vector<double> get_Linear_Transformation_expected_sequentialVector(int dimension, vector<vector<double>> input_matrix)
{
    vector<double> result(dimension);
    vector<double> startInput = input_matrix[0];
    vector<vector<double>> outputvector(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        if (i == 0)
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, startInput);
        }
        else
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, outputvector[i - 1]);
        }
    }
    int size = outputvector.size();
    return outputvector[size - 1];
}

inline void print_error_difference(vector<double> actualResult, vector<double> expectedResult, int dimention)
{
    for (unsigned int i = 0; i < dimention; i++)
    {
        cout << "index " << i << ": " << expectedResult[i] - actualResult[i] << "\n";
    }
}

inline void get_max_error_norm(vector<double> actualResult, vector<double> expectedResult, int dimension, vector<double> errordResult)
{
    for (int i = 0; i < dimension; i++)
    {
        errordResult[i] = abs(expectedResult[i] - actualResult[i]);
    }

    // Find the maximum element in error vector
    auto maxError = max_element(errordResult.begin(), errordResult.end());
    double indx = distance(errordResult.begin(), maxError);
    cout << "\nMax Error Norm = " << *maxError<< " At Index: "<<indx <<"\n";
}

//Normalised Max Error Norm
inline double get_max_error_norm_normalised(int highestMatrixValue, vector<double> actualResult, vector<double> expectedResult, int dimension, vector<double> errordResult)
{
    for (int i = 0; i < dimension; i++)
    {
        errordResult[i] = abs(expectedResult[i] - actualResult[i]);
    }
    cout << "\nMax Error Norm = " << *max_element(errordResult.begin(), errordResult.end()) << "\n";
    double maxNorm = *max_element(errordResult.begin(), errordResult.end()) / highestMatrixValue;
    return maxNorm;
}

inline double calc_mse(int dimension, vector<double> errordResult)
{
    double mseValue = 0.0;
    for (int j = 0; j < dimension; j++)
    {
        mseValue += errordResult[j] * errordResult[j];
    }
    return mseValue;
}

inline double get_mse(vector<double> actualResult, vector<double> expectedResult, int dimension, vector<double> errordResult)
{
    for (int i = 0; i < dimension; i++)
    {
        errordResult[i] = expectedResult[i] - actualResult[i];
    }
    double mseValue = calc_mse(dimension, errordResult);
    return mseValue;
}

//Test Phase
inline void get_max_error_norm_largeValues(vector<double> actualResult, vector<double> expectedResult, int dimension)
{
    vector<float> error(dimension);
    for (unsigned int i = 0; i < dimension; i++)
    {
        float exp = round(expectedResult[i] * 1000.0) / 1000.0;
        float act = round(actualResult[i] * 1000.0) / 1000.0;
        float difference = exp - act;
        error[i] = abs(difference);
    }

    // Find the maximum element in error vector
    cout << "\nMax Error Norm = " << *max_element(error.begin(), error.end()) << "\n";

}

inline void print_Linear_Transformation_expected_sequentialVector(int dimension, vector<vector<double>> input_matrix)
{
    vector<double> result(dimension);
    vector<double> startInput = input_matrix[0];
    vector<vector<double>> outputvector(dimension);
    int k = 0;
    for (int i = 0; i < dimension; i++)
    {
        if (i == 0)
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, startInput);
        }
        else
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, outputvector[i - 1]);
        }
    }
    int size = outputvector.size();
    /*vector<double> result = outputvector[size - 1];*/
    for (unsigned int row = 0; row < dimension; row++)
    {
        cout << "[";
        for (unsigned int x = 0; x < dimension; x++)
        {
            cout << outputvector[row][x] << ", ";
        }
        cout << "]" << endl;

    }
}

inline void print_Linear_Transformation_expected_sequential_dense(int iteration, int dimension, vector<vector<double>> input_matrix)
{
    cout << "\n    + Printing expected values:" << endl;
    vector<double> result(dimension);
    vector<double> startInput = input_matrix[0];
    vector<vector<double>> outputvector(dimension);
    int k = 0;
    for (int i = 0; i < iteration; i++)
    {
        if (i == 0)
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, startInput);
        }
        else
        {
            outputvector[i] = get_Linear_Transformation_expected_vector(dimension, input_matrix, outputvector[i - 1]);
        }
    }
    int size = outputvector.size();
    /*vector<double> result = outputvector[size - 1];*/
    for (unsigned int row = 0; row < iteration; row++)
    {
        cout << "\n    + Iteration:" << row << endl;
        cout << "[";
        for (unsigned int x = 0; x < dimension; x++)
        {
            cout << outputvector[row][x] << ", ";
        }
        cout << "]" << endl;

    }
}

inline Ciphertext Linear_Transform_Cipher_Sequential_2(Ciphertext ct, vector<Ciphertext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params,
    bool rescale, SecretKey sk, RelinKeys relin_keys, vector<Plaintext> plainDiagonal, PublicKey pk, vector<vector<double>> all_diagonal, double scale)
{
    cout << "    + Invoked Linear_Transform_Cipher C_vec . C_mat:" << endl;
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    Encryptor encryptor(context, pk);
    double scaleSet = pow(2.0, 80);
    // Create CKKS encoder
    CKKSEncoder ckks_encoder(context);

    int size = U_diagonals.size();
    vector<Ciphertext> cipher_result_prime(size);
    for (int i = 0; i < size; i++)
    {
        Ciphertext ct_rot;

        Ciphertext ct_new;
        if (i == 0)
        {
            evaluator.rotate_vector(ct, -size, gal_keys, ct_rot);
            evaluator.add(ct, ct_rot, ct_new);
        }
        else
        {
            evaluator.rotate_vector(cipher_result_prime[i - 1], -size, gal_keys, ct_rot);
            evaluator.add(cipher_result_prime[i - 1], ct_rot, ct_new);
        }

        cout << "    + Scale on add - input vector: " << log2(ct_new.scale()) << " bits" << endl;
        cout << "    + Scale on add - input matrix: " << log2(U_diagonals[0].scale()) << " bits" << endl;

        int ctNewCoeff = ct_new.coeff_modulus_size();
        int inputMatrixCoeff = U_diagonals[0].coeff_modulus_size();
        int differenceInCoeff = inputMatrixCoeff - ctNewCoeff;

        for (int p = 0; p < differenceInCoeff; p++)
        {
            cout << "    Mod switch to next inplace on InputMatrix : " << p << endl;
            for (int r = 0; r < size; r++)
            {
                evaluator.mod_switch_to_next_inplace(U_diagonals[r]);
            }
        }

        vector<Ciphertext> ct_result(size);
        evaluator.multiply(ct_new, U_diagonals[0], ct_result[0]);
        cout << "    + Scale on multiply - input vector : " << log2(ct_result[0].scale()) << " bits" << endl;
        cout << "    + Scale on multiply - input matrix: " << log2(U_diagonals[0].scale()) << " bits" << endl;

        /* cout << "    + Scale on multiply - input vector - before relin (1st element): " << log2(ct_result[0].scale()) << " bits" << endl;
         evaluator.relinearize_inplace(ct_result[0], relin_keys);
         cout << "    + Scale on multiply - input vector - after relin (1st element): " << log2(ct_result[0].scale()) << " bits" << endl;
         cout << "    + Scale on multiply - input vector - before rescale (1st element): " << log2(ct_result[0].scale()) << " bits" << endl;
         evaluator.rescale_to_next_inplace(ct_result[0]);
         cout << "    + Scale on multiply - input vector - after rescale (1st element): " << log2(ct_result[0].scale()) << " bits" << endl;

         cout << "    + Scale on multiply - input matrix: " << log2(U_diagonals[0].scale()) << " bits" << endl;*/
        for (int l = 1; l < size; l++)
        {
            Ciphertext temp_rot;
            evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
            evaluator.multiply(temp_rot, U_diagonals[l], ct_result[l]);
        }

        //for (int j = 1; j < size; j++)
        //{            
        //    // Added relin and rescale
        //    evaluator.relinearize_inplace(ct_result[j], relin_keys);
        //    evaluator.rescale_to_next_inplace(ct_result[j]);
        //}


        cout << "    + Scale on multiply - input matrix - after and relin rescale (all elements): " << log2(ct_result[1].scale()) << " bits" << endl;

        evaluator.add_many(ct_result, cipher_result_prime[i]);

        cout << "    + Scale on add_many - input vector - before relin: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
        evaluator.relinearize_inplace(cipher_result_prime[i], relin_keys);
        cout << "    + Scale on add_many - input vector - after relin: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;

        if (rescale)
        {
            cout << "    + Scale on add_many - input vector - before rescale: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
            evaluator.rescale_to_next_inplace(cipher_result_prime[i]);
            cout << "    + Scale on add_many - input vector - after rescale: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
        }

        Plaintext pt_result_current;
        decryptor.decrypt(cipher_result_prime[i], pt_result_current);

        // Decode
        vector<double> output_result_plain;
        ckks_encoder.decode(pt_result_current, output_result_plain);

        cout << "Linear Transformation Set 2 Result:" << endl;
        /*print_partial_vector(output_result_plain, size); */
        for (unsigned int row = 0; row < size; row++)
        {
            cout << output_result_plain[row] << ", ";
        }
    }

    return cipher_result_prime[size - 1];
}

inline Ciphertext Linear_Transform_Cipher_Sequential_Dense(Ciphertext ct, vector<Ciphertext> U_diagonals, GaloisKeys gal_keys, EncryptionParameters params,
    bool rescale, SecretKey sk, RelinKeys relin_keys, int iteration)
{
    cout << "    + Invoked Linear_Transform_Cipher C_vec . C_mat:" << endl;
    auto context = SEALContext::SEALContext(params);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    CKKSEncoder ckks_encoder(context);
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    int size = U_diagonals.size();
    vector<Ciphertext> cipher_result_prime(size);
    for (int i = 0; i < iteration; i++)
    {
        cout << "\nIteration: " << i << endl;

        /*Initialization*/
        auto start = chrono::high_resolution_clock::now();
        Ciphertext ct_rot;
        Ciphertext ct_new;

        /*Rotate and Add*/
        if (i == 0)
        {
            evaluator.rotate_vector(ct, -size, gal_keys, ct_rot);
            evaluator.add(ct, ct_rot, ct_new);
        }
        else
        {
            evaluator.rotate_vector(cipher_result_prime[i - 1], -size, gal_keys, ct_rot);
            evaluator.add(cipher_result_prime[i - 1], ct_rot, ct_new);
        }

        cout << "\n    + Scale on add - input vector: " << log2(ct_new.scale()) << " bits" << endl;
        cout << "\n    + Scale on add - input matrix: " << log2(U_diagonals[0].scale()) << " bits" << endl;

        /*Mod inplace on input matrix*/
        int ctNewCoeff = ct_new.coeff_modulus_size();
        int inputMatrixCoeff = U_diagonals[0].coeff_modulus_size();
        int differenceInCoeff = inputMatrixCoeff - ctNewCoeff;
        for (int p = 0; p < differenceInCoeff; p++)
        {
            cout << "    Mod switch to next inplace on InputMatrix : " << p << endl;
            for (int r = 0; r < size; r++)
            {
                evaluator.mod_switch_to_next_inplace(U_diagonals[r]);
            }
        }

        /*Multiply with all rows*/
        vector<Ciphertext> ct_result(size);
        evaluator.multiply(ct_new, U_diagonals[0], ct_result[0]);
        cout << "\n    + Scale on multiply - input vector : " << log2(ct_result[0].scale()) << " bits" << endl;
        cout << "\n  + Scale on multiply - input matrix: " << log2(U_diagonals[0].scale()) << " bits" << endl;
        for (int l = 1; l < size; l++)
        {
            Ciphertext temp_rot;
            evaluator.rotate_vector(ct_new, l, gal_keys, temp_rot);
            evaluator.multiply(temp_rot, U_diagonals[l], ct_result[l]);
        }

        /*Add many, relin, rescale*/
        evaluator.add_many(ct_result, cipher_result_prime[i]);
        cout << "\n    + Scale on add_many - input vector - before relin: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
        evaluator.relinearize_inplace(cipher_result_prime[i], relin_keys);
        cout << "\n    + Scale on add_many - input vector - after relin: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
        if (rescale)
        {
            cout << "\n    + Scale on add_many - input vector - before rescale: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
            evaluator.rescale_to_next_inplace(cipher_result_prime[i]);
            cout << "\n    + Scale on add_many - input vector - after rescale: " << log2(cipher_result_prime[i].scale()) << " bits" << endl;
        }
        auto stop = chrono::high_resolution_clock::now();

        /*Print time*/
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "\nTime to compute Iteration: " << i << duration.count() << " microseconds" << endl;


        cout << "\nVirtual Memory to compute Iteration: " << i + 1 << pmc.PrivateUsage << endl;
        cout << "\nPhysical Memory (RAM) to compute Iteration: " << i + 1 << pmc.WorkingSetSize << endl;

        /*Test with decrypt - to be removed*/
        Plaintext pt_result_current;
        decryptor.decrypt(cipher_result_prime[i], pt_result_current);
        vector<double> output_result_plain;
        ckks_encoder.decode(pt_result_current, output_result_plain);
        cout << "\nLinear Transformation Set 2 Result:" << endl;
        for (unsigned int row = 0; row < size; row++)
        {
            cout << output_result_plain[row] << ", ";
        }
    }

    return cipher_result_prime[size - 1];
}



//Filewriter

inline void SaveInFile(string completePath, string content)
{
    std::ofstream myfile;
    myfile.open(completePath, std::ios::app);
    myfile << content;
    myfile.close();
}