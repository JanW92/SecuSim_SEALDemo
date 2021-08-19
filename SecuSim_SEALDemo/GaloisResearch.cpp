#include "CoreFunctionality.h"


void InvestigateGaloisKeySize()
{
    EncryptionParameters parms(scheme_type::ckks);
    
    //size_t poly_modulus_degree = 8192;
    size_t poly_modulus_degree = 32768;

    /*
        Precision of 10 ^ -7 requires at least scale of 23 bit but noise must be concidered as well.
        Resolution of the number before the dot must be concidered as well: 8 < N -> 3 bit
        first prime - scale = bits for number before dot
        60 -40 = 20
    */
    double scale = pow(2.0, 40);
    /*
    +----------------------------------------------------+
    |poly_modulus_degree | max coeff_modulus bit - length|
    +---------------------+------------------------------+
    | 1024 | 27                                             |
    | 2048 | 54                                             |
    | 4096 | 109                          |
    | 8192 | 218                          |
    | 16384 | 438                          |
    | 32768 | 881 |
    +-------------------- - +------------------------------ +
    */
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60 }));
    
    SEALContext context(parms);
    print_parameters(context);
    
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    
    const int rotations = 10;
    vector<int> rotationVect{ 1 };

    GaloisKeys gal_keys;
    keygen.create_galois_keys(rotationVect,gal_keys);
    //keygen.create_galois_keys( gal_keys);

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    CKKSEncoder encoder(context);

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
        //expectedSolution[i] = toBeRotated[i] * scalingVector[i];
        expectedSolution[i] = toBeRotated[i];
    }
    

    Plaintext plainScale, plaineRot, plainSolution;
    encoder.encode(toBeRotated,scale, plaineRot);
    encoder.encode(scalingVector, scale, plainScale);

    Ciphertext cipherScale, cipcherRot,cipherSolution;
    encryptor.encrypt(plainScale, cipherScale);
    encryptor.encrypt(plaineRot, cipcherRot);

    //evaluator.multiply (cipherScale, cipcherRot, cipherSolution);
    
    for (int i = 0; i < rotations +1; i++)
    {
        //Rotate encrypted and unencrypted vector
        rotate(expectedSolution.begin(), expectedSolution.begin() + 1, expectedSolution.end());
        evaluator.rotate_vector_inplace(cipcherRot, 1, gal_keys);

        decryptor.decrypt(cipcherRot, plainSolution);
        vector<double> solution;
        encoder.decode(plainSolution, solution);
        std::cout <<  i << ". rotation:\n";
        get_max_error_norm(solution, expectedSolution, expectedSolution.size(),errorResult);
    }


}