main:
	g++ -std=c++14 main.cpp
	./a.out

git:
	git add .
	git commit -m "$m"
	git push -u origin master