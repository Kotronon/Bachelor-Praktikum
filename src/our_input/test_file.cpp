//
// Created by Samsamsamsamsamsamsamsamsa on 03.12.2023.
//

#include "test_file.h"
//#include "newinput-pimpl.cxx"

#include <string>
#include <iostream>

//alle parameter!


int test;


void getint(int a){
        test = a;
}



int main(int argc, char *argsv[]){

    std::cout<<"Enter number"<< std::endl;
    std::cin>>test;

    getint(test);

    std::cout<<"Your name was "<<test<< std::endl;


}
