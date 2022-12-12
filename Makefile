main:
	g++ -std=c++17 *.cpp
	./a.out

git:
	git add .
	git commit -m "$m"
	git push -u origin master