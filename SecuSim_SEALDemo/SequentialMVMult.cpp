#include "CoreFunctionality.h"
#include <chrono>
#include <omp.h>

void CreateMatrix(vector<vector<double>>& matrix);
vector<double>  MatrixVectorMultiplication(vector<vector<double>> matrix, vector<double> fieldVector);


void InvestigateMVPerformance()
{
    std::cout << std::scientific;
    std::cout << "Begin of encrypted Matrix-Vektor multiplication\n";
    //omp_set_num_threads(4);
    string fileName = "SecuSim_SEALMeasurement.txt";
    string path = "D:\\PerformanceAnalyse\\";
    EncryptionParameters parms(scheme_type::ckks);

    /*
    *   Maximum  poly_modulus_degree
    */
    size_t poly_modulus_degree = 8192;
    /*
        The matrix vector multiplication with poly_modulus_degree 32768 or 16384 is not feasible. For the corresponding matrix are 2 respectively 0,5 GB of memory required. 
        Encryption comes with an overhead of ~150 respectively ~50 for the maximum poly_modulus_degree and corresponding max coeff_modulus (see table in project description).
        This results in 300 GB respectively 25 GB memory requirement.
    */
    //size_t poly_modulus_degree = 32768;
    /*
    * The scale value defines the resolution of the digits after the dot (the float number is multiplied with the scale number to get an integer that is used furhter on)
    *  ->  Precision of 10 ^ -7 requires at least scale of 23 bit but noise must be concidered as well. The exactly required amount of bits to cover the noise is not possible (really hard) to calculate.
    *  ->  Resolution of the number before the dot must be concidered as well: Max value = 8 -> N mind 3 bit
    *  -> first prime - scale = bits for number before dot: 60 - 40 = 20
    *  -> Remaining primes are used for rescaling and schould have the same size as the scale value
    */

    double scale = pow(2.0, 40);
    /*
    * Recommanded poly_modulus_degrees and corresponding max coeff_modulus bit length
    *
    *
    *     +-----------------------------------------------------+
    *     |poly_modulus_degree  | max coeff_modulus bit - length|
    *     +---------------------+-------------------------------+
    *     | 1024                | 27                            |
    *     | 2048                | 54                            |
    *     | 4096                | 109                           |
    *     | 8192                | 218                           |
    *     | 16384               | 438                           |
    *     | 32768               | 881                           |
    *     +-------------------- - +---------------------------- +
    */

    parms.set_poly_modulus_degree(poly_modulus_degree);

    /*
        Amount of primes between the firstand last are limiting the amount of operations!In this case 19 operations are possible.
    */
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    //parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60 }));
    //parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 60 }));
    //const int operations = 18;
    const int operations = 2;

    /*
        Setup SEAL objects
    */
    SEALContext context(parms);
    print_parameters(context);
    string contextAsString = get_parameters(context);
    SaveInFile(path + fileName, contextAsString);


    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    CKKSEncoder encoder(context);

    /*
        Reduced galois key that supports rotations by one in the leftern direction.For details see function : InvestigateGaloisKeySize
     */
    vector<int> rotationVect{ 1 };
    GaloisKeys gal_keys;
    keygen.create_galois_keys(rotationVect, gal_keys);

    /*
        Creation of vectors.Their size is predertermined by the poly_modulus_degree value.Since the ckks scheme is usedand it uses only real numbers the vector size is calculated as followed :
    */
    size_t slot_count = encoder.slot_count();

    vector<double> temp_solution;
    vector<double> fieldDummy(slot_count), solution(slot_count);
    vector<double> errorResult(slot_count);

    auto start_matrixCreation = chrono::high_resolution_clock::now();

    /*
    * Creation Matrix
    */
    vector<vector<double>> curlDummy(slot_count, vector<double>(slot_count));
    //vector<vector<double>> matrixDia(slot_count, vector<double>(slot_count));
    vector<vector<double>> matrixDia(slot_count, vector<double>(slot_count));

    CreateMatrix(curlDummy);

    //#pragma omp parallel
    //{
    //    #pragma omp for
        for (int i = 0; i < slot_count; i++)
        {
            matrixDia[i] = get_diagonal(i, curlDummy);
        }
 //   }
    curlDummy.clear();
    curlDummy.shrink_to_fit();

    auto after_matrixCreation = chrono::high_resolution_clock::now();

    /*
        Vector initialization
    */

    double value = 2.01234567890123456789;

    for (int i = 0; i < slot_count; i++)
    {
        fieldDummy[i] = value;
        solution[i] = 0.;
    }

    /*
        Calculate expected result after one matrix-vector multiplication
    */
    auto before_PlainMV = chrono::high_resolution_clock::now();
    vector<double> expectedSolution = MatrixVectorMultiplication(matrixDia, fieldDummy);
    auto after_PlainMV = chrono::high_resolution_clock::now();

    /*
        Canonial embedding of the data vectors
    */
    Plaintext p_Field, p_Temp, p_Solution;
    encoder.encode(fieldDummy, scale, p_Field);
    encoder.encode(solution, scale, p_Solution);

    /*
        Encryption
    */
    auto before_Encryption = chrono::high_resolution_clock::now();
    vector<Ciphertext> c_Matrix(slot_count);
    Ciphertext c_Field, c_Solution, c_NullVector;
    encryptor.encrypt(p_Field, c_Field);
    encryptor.encrypt(p_Solution, c_Solution);
    /*Result of multiplication is going to be added with c_Solution. Therefore the mod has to be adjusted.*/
    evaluator.mod_switch_to_next_inplace(c_Solution);
    c_NullVector = c_Solution;

    //#pragma omp parallel
    //{
    //    #pragma omp for
        for (int i = 0; i < slot_count; i++)
        {
            encoder.encode(matrixDia[i], scale, p_Temp);
            encryptor.encrypt(p_Temp, c_Matrix[i]);
            p_Temp.release();
        }
//    }

    //decryptor.decrypt(c_Matrix[0], p_Temp);
    //encoder.decode(p_Temp, solution);

    //decryptor.decrypt(c_Matrix[1], p_Temp);
    //encoder.decode(p_Temp, solution);

    //decryptor.decrypt(c_Field, p_Temp);
    //encoder.decode(p_Temp, solution);

    Ciphertext c_Temp;
    auto after_Encryption = chrono::high_resolution_clock::now();
    vector<double> operationTimings(operations);
    vector<double> operationPrecision(operations);
    /*
        Matrix-Vector multiplikation
    */

    //#pragma omp parallel
    //{
    //  #pragma omp for
        for (int j = 1; j <= operations; j++)
        {
            auto before_EncryptedMV = chrono::high_resolution_clock::now();
            for (int i = 0; i < c_Matrix.size(); i++)
            {

                    evaluator.multiply(c_Matrix[i], c_Field, c_Temp);
                    evaluator.relinearize_inplace(c_Temp, relin_keys);
                    evaluator.rescale_to_next_inplace(c_Temp);

                    //decryptor.decrypt(c_Temp, p_Temp);
                    //encoder.decode(p_Temp, temp_solution);

                    c_Solution.scale() = scale;
                    c_Temp.scale() = scale;
                    evaluator.add_inplace(c_Solution, c_Temp);
                    c_Temp.release();

                    //decryptor.decrypt(c_Solution, p_Temp);
                    //encoder.decode(p_Temp, temp_solution);
                    evaluator.rotate_vector_inplace(c_Field, 1, gal_keys);

                    evaluator.mod_switch_to_next_inplace(c_Matrix[i]);
            }
            auto after_EncryptMV = chrono::high_resolution_clock::now();
            /*
                Compare the results
            */
            c_Field = c_Solution;
            decryptor.decrypt(c_Field, p_Temp);
            encoder.decode(p_Temp, temp_solution);
            operationPrecision[j-1] = get_max_error_norm_value(temp_solution, expectedSolution, expectedSolution.size(), errorResult);
            auto t4 = std::chrono::duration_cast<std::chrono::duration<double>>(after_EncryptMV - before_EncryptedMV);
            operationTimings[j-1] = t4.count();
            /*
                Preperations for the next operation: Update expected result and set c_Solution back to a zero filled vector with reduced mod
            */
            if (j < operations)
            {
                expectedSolution = MatrixVectorMultiplication(matrixDia, expectedSolution);
                evaluator.mod_switch_to_next_inplace(c_NullVector);
                c_Solution = c_NullVector;
            }
        }

        //auto t1 = std::chrono::duration_cast<std::chrono::microseconds>(after_matrixCreation - start_matrixCreation);
        auto t1 = std::chrono::duration_cast<std::chrono::duration<double>>(after_matrixCreation - start_matrixCreation);
        auto t2 = std::chrono::duration_cast<std::chrono::duration<double>>(after_PlainMV - before_PlainMV);
        auto t3 = std::chrono::duration_cast<std::chrono::duration<double>>(after_Encryption - before_Encryption);
        double totalTime(0);

        std::ostringstream streamObj;
        
        string output= "Time for matrix assembling: " + to_string(t1.count()) + " s\n"
            + "Time for encryption: " + to_string(t3.count()) + " s\n";
        for (int i = 0; i < operations; i++)
        {
            streamObj << operationPrecision[i];

            totalTime += operationTimings[i];
            output += to_string(i) + ". operation:\n" +
                + "Time for encrypted MV: " + to_string(operationTimings[i]) + " s\n"
                + "Lowest precision: " + streamObj.str() + "\n";
            streamObj.str("");

        }

        /*
            An EMW simulation with the electric (E) and magnetic (H) field components as degrees of freedom in a 3 dimensional space contains 6 values per node (Ex,Ey,Ez,Hx,Hy,Hz).
            Thus a matrix with dimension nxn calculates n/6 nodes.
        */
        int nodes = std::floor(slot_count / 6.);

        /*
            Performance measured in Mega Cells per Second [MCells/s]
        */
        double performance = (nodes * operations) / (totalTime * pow(10, 6));
        output += "Measured simulation performance: " + to_string(performance) + " MCells/sec";
        std::cout << output;
        SaveInFile(path + fileName, output);
  //  }
}

/**********************************************************************/
/**************************** Helper function *************************/
/**********************************************************************/

vector<double>  MatrixVectorMultiplication(vector<vector<double>> matrix, vector<double> fieldVector)
{
    vector<double> resultVector(matrix.size(),0.);

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
/*
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
*/

/*
* Creates a scale matrix.
*/
void CreateMatrix(vector<vector<double>>& matrix)
{

    //double val = 0.55;
    //for (int i = 0; i < matrix.size(); i++)
    //{
    //    for (int j = 0; j < matrix.size(); j++)
    //    {
    //        if (j == i)
    //        {
    //            matrix[i][j] = val;
    //        }
    //        else
    //        {
    //            matrix[i][j] = 0.;
    //        }
    //    }
    //}
    double val = 1.1;
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
                matrix[i][j] = 0.;
            }
        }
    }
}
