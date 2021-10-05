#include "CoreFunctionality.h"

void ScaleVector(vector<double>& values, double scale);

void InvestigateGaloisKeySize()
{
    EncryptionParameters parms(scheme_type::ckks);
    
    //size_t poly_modulus_degree = 8192;
    size_t poly_modulus_degree = 32768;


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
    //The CoeffModulus has a huge impact on the size of the SEALContext object, which is the largest object of them all. In the maximum case with poly_modulus_degree = 32768 the object is ~500 MB big
    //parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60 }));
    

    SEALContext context(parms);
    print_parameters(context);
    
    //SEAL objects
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
    

    //The amount of rotations seeems not to be limited, since it does not require rescaling nor relinearization
    const int rotations = 10;

    //The amount of steps and the direction of the  rotation are determined by the galois key. 
    //The default constructor makes rotations in all directions and by many different step sizes possible. This as well as the poly_modulus_degree have a big impact on the galois key. 
    //The maximum case generates a galois key with 6 GB !
    //keygen.create_galois_keys( gal_keys);
    
    //The alternative is a reduced key that can only perform rotations in a given step size and direction.
    //By passing rotationVect in the constructor the galois key can perform rotations to the left by 1
    vector<int> rotationVect{ 1 };
    GaloisKeys gal_keys;
    keygen.create_galois_keys(rotationVect,gal_keys);


    vector<double> toBeRotated, scalingVector, expectedSolution, errorResult;
    size_t slot_count = encoder.slot_count();
    double value = 1.01234567890123456789;
    double scaleVal = 1.5;
    toBeRotated.resize(slot_count);
    scalingVector.resize(slot_count);
    expectedSolution.resize(slot_count);
    errorResult.resize(slot_count);

    for (int i = 0; i < slot_count; i++)
    {
        scalingVector[i] = scaleVal;
        toBeRotated[i] = value + i%4;
        expectedSolution[i] = toBeRotated[i];
    }
    

    Plaintext plainScale, plaineRot, plainSolution;
    encoder.encode(toBeRotated,scale, plaineRot);
    encoder.encode(scalingVector, scale, plainScale);

    Ciphertext cipherScale, cipcherRot,cipherSolution;
    encryptor.encrypt(plainScale, cipherScale);
    encryptor.encrypt(plaineRot, cipcherRot);
    
    for (int i = 0; i < rotations +1; i++)
    {
        //Rotate encrypted and unencrypted vector
        rotate(expectedSolution.begin(), expectedSolution.begin() + 1, expectedSolution.end());
        evaluator.rotate_vector_inplace(cipcherRot, 1, gal_keys);

        //Compare results
        decryptor.decrypt(cipcherRot, plainSolution);
        vector<double> solution;
        encoder.decode(plainSolution, solution);
        std::cout <<  i << ". rotation:\n";
        get_max_error_norm(solution, expectedSolution, expectedSolution.size(),errorResult);
    }


}

void InvestigateRotationSequentialOperation()
{
    EncryptionParameters parms(scheme_type::ckks);
    
    //Maximum  poly_modulus_degree
    size_t poly_modulus_degree = 32768;

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
    
    //Amount of primes between the first and last are limiting the amount of operations! In this case 19 operations are possible.
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60 }));
    const int operations = 19;

    //Setup SEAL objects
    SEALContext context(parms);
    print_parameters(context);

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

    //Reduced galois key that supports rotations by one in the leftern direction. For details see function: InvestigateGaloisKeySize
    vector<int> rotationVect{ 1 };
    GaloisKeys gal_keys;
    keygen.create_galois_keys(rotationVect, gal_keys);

    //Creation of vectors. Their size is predertermined by the poly_modulus_degree value. Since the ckks scheme is used and it uses only real numbers the vector size is calculated as followed:
    //slot_count = poly_modulus_degree / 2
    vector<double> toBeRotated, scalingVector, expectedSolution, errorResult;
    size_t slot_count = encoder.slot_count();
    toBeRotated.resize(slot_count);
    scalingVector.resize(slot_count);
    expectedSolution.resize(slot_count);
    errorResult.resize(slot_count);

    //Pay attention that the final value doesn't grow too much. If scaleVal is set to 1.5 the accuracy drops after 8 itterations, 
    //since the 20 bits of pre dot resolution is not sufficient enough anymore.
    double scaleVal = 1.1;
    double value = 1.01234567890123456789;

    for (int i = 0; i < slot_count; i++)
    {
        scalingVector[i] = scaleVal;
        toBeRotated[i] = value + i % 4;
        expectedSolution[i] = toBeRotated[i];
    }

    Plaintext plainScale, plaineRot, plainSolution;
    encoder.encode(toBeRotated, scale, plaineRot);
    encoder.encode(scalingVector, scale, plainScale);

    Ciphertext cipherScale, cipherRot, cipherSolution, tempSol;
    encryptor.encrypt(plainScale, cipherScale);
    encryptor.encrypt(plaineRot, cipherRot);


    for (int i = 0; i < operations ; i++)
    {

        //Rotate and scale unencrypted vector
        rotate(expectedSolution.begin(), expectedSolution.begin() + 1, expectedSolution.end());
        ScaleVector(expectedSolution, scaleVal);
        
        //Rotate and scale encrypted vector
        evaluator.rotate_vector_inplace(cipherRot, 1, gal_keys);
        evaluator.multiply_inplace(cipherRot, cipherScale);

        //A multiplication increases the size of a ciphertext. By relinearization it is reduced to its original size.
        evaluator.relinearize_inplace(cipherRot, relin_keys);
        //The scale squares in a multiplication. Therefor the ciphertext has to be rescaled, which consumes one prime of the CoeffModulus
        evaluator.rescale_to_next_inplace(cipherRot); 
        //Only vectors with the same amount of primes can be multiplied. Thus the cipherScale needs to loose one prime, without reducing the scale. Rescaling would decrese the scale to 1.
        //For such a case mod_switch can be used to select the next active prime
        evaluator.mod_switch_to_next_inplace(cipherScale);

        //Compare the results
        decryptor.decrypt(cipherRot, plainSolution);
        vector<double> solution;
        encoder.decode(plainSolution, solution);
        std::cout << i << ". rotation:\n";
        get_max_error_norm(solution, expectedSolution, expectedSolution.size(), errorResult);
    }

}

/**********************************************************************/
/**************************** Helper function *************************/
/**********************************************************************/
void ScaleVector(vector<double> &values, double scale)
{
    for (int i = 0; i < values.size(); i++)
    {
        values[i] *= scale;
    }
}