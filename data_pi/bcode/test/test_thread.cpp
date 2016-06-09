#include <iostream>
#include <thread>

void call_thread(){
	std::cout << "something" << std::endl;
}

int main(){
	std::thread t1(call_thread);
	t1.join();
	return 0;
}
