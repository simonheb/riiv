#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;
#include <tictoc.h>


mat tictoc_timers=zeros(100,2);
int tictoc_current=0;
auto start=std::chrono::system_clock::now();

void toc() {
  //  std::cout<<"ending" << current <<endl;
  tictoc_timers(tictoc_current,1)+=std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::system_clock::now() - start).count();
  tictoc_current = 0;
  tictoc_timers.col(0)=linspace(0,99);
}
void tic(int counter){ 
  toc();
  tictoc_current = counter;
  //std::cout<<"starting" << current <<endl;
  start=std::chrono::system_clock::now();
}
void tictoc(int limit){ 
  toc();
  std::cout<<tictoc_timers.head_rows(limit) <<endl;
}