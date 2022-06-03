#include <iostream>
#include <unistd.h>
#include "Point.h"
#include "Strategy.h"
#include "Judge.h"
#include <vector>
// #include <sys/time.h>
#include <cmath>
#include <ctime>

using namespace std;

/*
	策略函数接口,该函数被对抗平台调用,每次传入当前状态,要求输出你的落子点,该落子点必须是一个符合游戏规则的落子点,不然对抗平台会直接认为你的程序有误
	
	input:
		为了防止对对抗平台维护的数据造成更改，所有传入的参数均为const属性
		M, N : 棋盘大小 M - 行数 N - 列数 均从0开始计， 左上角为坐标原点，行用x标记，列用y标记
		top : 当前棋盘每一列列顶的实际位置. e.g. 第i列为空,则_top[i] == M, 第i列已满,则_top[i] == 0
		_board : 棋盘的一维数组表示, 为了方便使用，在该函数刚开始处，我们已经将其转化为了二维数组board
				你只需直接使用board即可，左上角为坐标原点，数组从[0][0]开始计(不是[1][1])
				board[x][y]表示第x行、第y列的点(从0开始计)
				board[x][y] == 0/1/2 分别对应(x,y)处 无落子/有用户的子/有程序的子,不可落子点处的值也为0
		lastX, lastY : 对方上一次落子的位置, 你可能不需要该参数，也可能需要的不仅仅是对方一步的
				落子位置，这时你可以在自己的程序中记录对方连续多步的落子位置，这完全取决于你自己的策略
		noX, noY : 棋盘上的不可落子点(注:涫嫡饫锔?龅膖op已经替你处理了不可落子点，也就是说如果某一步
				所落的子的上面恰是不可落子点，那么UI工程中的代码就已经将该列的top值又进行了一次减一操作，
				所以在你的代码中也可以根本不使用noX和noY这两个参数，完全认为top数组就是当前每列的顶部即可,
				当然如果你想使用lastX,lastY参数，有可能就要同时考虑noX和noY了)
		以上参数实际上包含了当前状态(M N _top _board)以及历史信息(lastX lastY),你要做的就是在这些信息下给出尽可能明智的落子点
	output:
		你的落子点Point
*/

// #define TIME 2500000
#define TIME 2.5

class Naive{
private:
	int M, N, lastX, lastY, noX, noY;
	const int *top;
	int **board;
public:
	Naive(const int M, const int N, const int *top, int **board, const int lastX, const int lastY, const int noX, const int noY): M(M), N(N), top(top), board(board), lastX(lastX), lastY(lastY), noX(noX), noY(noY){}
	~Naive(){}
	Point getPoint();
};

struct Node{
	int player, winner; // 0,1,2分别对应无落子,有用户的子,有程序的子
	int Nodewins, T;
	double I, X; // UCB1的信心上界索引
	std::vector<int> top;
	std::vector<Node*> children;
	Node(int player, int N, int M, int **board): player(player), winner(0), Nodewins(0), T(0), I(0.0), X(0.0){
        children.resize(N);
        top.resize(N);
        for(int i = 0; i < N; i++){
            int x = M - 1;
            while(x >= 0 && board[x][i] != 0) x--;
            children[i] = nullptr;
            top[i] = x;
        }
    }
	~Node(){
		for(int i = 0; i < children.size(); i++){
			if(children[i]) delete children[i];
		}
	}
};

class UCT{
private:
	int M, N, lastX, lastY, noX, noY;
	const int *top;
	int **board, **baseBoard;
	Node *root;
protected:
    Node* treePolicy(std::vector<Node*>& path);
	int simulation(Node *node);
    void backPropagation(std::vector<Node*>& path, int result);
	bool haveWin(int x, int y);

	int expand(Node* &current);
	int bestChild(Node* &current, Node* &next);
	Point bestChild();
	int defaultPolicy(Node *node);
	void backup(Node* &pos, std::vector<Node*>& path, int result);

public:
	UCT(const int M, const int N, const int *top, int **_board, const int lastX, const int lastY, const int noX, const int noY): M(M), N(N), top(top), lastX(lastX), lastY(lastY), noX(noX), noY(noY){
		board = new int *[M];
		baseBoard = new int *[M];
		for(int i = 0; i < M; i++){
			board[i] = new int[N];
			baseBoard[i] = new int[N];
			for(int j = 0 ; j < N; j++){
				board[i][j] = _board[i][j];
				baseBoard[i][j]= _board[i][j];
			}
		}
		board[noX][noY] = 3;
		baseBoard[noX][noY] = 3;
		
		root = new Node(2, N, M, board); // 2为machine
	};
	~UCT(){
		if(root) delete root;
		if(board) clearArray(M, N, board);
		if(baseBoard) clearArray(M, N, baseBoard);
	}
	Point getPoint();
};

// struct timeval startTime;
clock_t startTime;

extern "C" Point *getPoint(const int M, const int N, const int *top, const int *_board,
						   const int lastX, const int lastY, const int noX, const int noY)
{
    // gettimeofday(&startTime, NULL); // 计时开始
	startTime = clock(); // 计时开始

	/*
		不要更改这段代码
	*/
	int x = -1, y = -1; //最终将你的落子点存到x,y中
	int **board = new int *[M];
	for (int i = 0; i < M; i++)
	{
		board[i] = new int[N];
		for (int j = 0; j < N; j++)
		{
			board[i][j] = _board[i * N + j];
		}
	}

	/*
		根据你自己的策略来返回落子点,也就是根据你的策略完成对x,y的赋值
		该部分对参数使用没有限制，为了方便实现，你可以定义自己新的类、.h文件、.cpp文件
	*/
	//Add your own code below

	// Naive 
	Naive naive(M, N, top, board, lastX, lastY, noX, noY);
	Point naive_point = naive.getPoint();
	x = naive_point.x;
	y = naive_point.y;
	if(x != -1 && y != -1){
		clearArray(M, N, board);
		return new Point(x, y);
	}

	// UCT
	UCT uct(M, N, top, board, lastX, lastY, noX, noY);
	Point uct_point = uct.getPoint();
	x = uct_point.x;
	y = uct_point.y;

	/*
		不要更改这段代码
	*/
	clearArray(M, N, board);
	return new Point(x, y);
}

/*
	getPoint函数返回的Point指针是在本so模块中声明的，为避免产生堆错误，应在外部调用本so中的
	函数来释放空间，而不应该在外部直接delete
*/
extern "C" void clearPoint(Point *p)
{
	delete p;
	return;
}

/*
	清除top和board数组
*/
void clearArray(int M, int N, int **board)
{
	for (int i = 0; i < M; i++)
	{
		delete[] board[i];
	}
	delete[] board;
}

/*
	添加你自己的辅助函数，你可以声明自己的类、函数，添加新的.h .cpp文件来辅助实现你的想法
*/
Point Naive::getPoint(){
	int x = -1, y = -1;
	// 自己必赢
	for (int i = N - 1; i >= 0; i--) {
		if (top[i] > 0) {
			x = top[i] - 1;
			y = i;	
			// 简单策略
			board[x][y] = 2;
			if(machineWin(x, y, M, N, board)) return Point(x, y);
			board[x][y] = 0;
		}
	}
	
	// 对方必赢
	for (int i = N - 1; i >= 0; i--) {
		if (top[i] > 0) {
			x = top[i] - 1;
			y = i;
			// 简单策略
			board[x][y] = 1;
			if(userWin(x, y, M, N, board)) return Point(x, y);
			board[x][y] = 0;
		}
	}

	return Point(-1, -1);
}

bool UCT::haveWin(int x, int y){
	if(board[x][y] == 1) return userWin(x, y, M, N, board);
	return machineWin(x, y, M, N, board);
}

int UCT::expand(Node* &current){
	int expand_rank = -1;
	for(int i = 0; i < N; i++){
		if(current->top[i] >= 0 && current->children[i] == nullptr){
			expand_rank = i;
			break;
		}
	}
	return expand_rank;
}

int UCT::bestChild(Node* &current, Node* &next){
	double I = 0.0; // UCB1的信心上界索引
	int best_child_rank = -1;
	for(int i = 0; i < N; i++){
		if(current->children[i] != nullptr && I < current->children[i]->I){
			I = current->children[i]->I;
			next = current->children[i];
			best_child_rank = i;
		}
	}
	return best_child_rank;
}

Node* UCT::treePolicy(std::vector<Node*> &path){
    Node *current = root;
	path.push_back(current);
    while(current->winner == 0){
		for(int i = 0; i < N; i++){ // 启发式优化：有必胜策略就停止模拟
            if(current->children[i] != nullptr && current->children[i]->winner == current->player){
                current = current->children[i];
                current->T += 2;
                current->Nodewins += 2;
				return current;
            }
        }
		int expand_rank = expand(current);
		/*
        int expand_rank = -1;
        for(int i = 0; i < N; i++){
            if(current->top[i] >= 0 && current->children[i] == nullptr){
                expand_rank = i;
                break;
            }
		}
		*/
		if(expand_rank != -1){
			int x = current->top[expand_rank], y = expand_rank;
            board[x][y] = current->player;
            int next_player = (current->player == 1 ? 2 : 1); // 换player
            current->children[expand_rank] = new Node(next_player, N, M, board);
            current->children[expand_rank]->T = 2;
            if(haveWin(x, y)){
                current->children[expand_rank]->winner = current->player;
                current->children[expand_rank]->Nodewins = 2 * (current->player == 2 ? 1 : 0);
            }
            else current = current->children[expand_rank];
			return current;
		}
		else{	
            Node *next = nullptr;
			int best_child_rank = bestChild(current, next);
            if(next == nullptr) return next;
            else{
                int x = current->top[best_child_rank], y = best_child_rank;
                board[x][y] = current->player;
                current = next;
            }
		}
		/*
        if(expand_rank != -1){
            int x = current->top[expand_rank];
            int y = expand_rank;
            board[x][y] = current->player;
            int next_player = (current->player == 1 ? 2 : 1); // 换player
            current->children[expand_rank] = new Node(next_player, N, M, board);
            current->children[expand_rank]->T = 2;
            if(haveWin(x, y)){
                current->children[expand_rank]->winner = current->player;
                current->children[expand_rank]->Nodewins = 2 * (current->player == 2 ? 1 : 0);
            }
            else current = current->children[expand_rank];
			break;
        }
        else{
            double I = 0.0; // UCB1的信心上界索引
            int road_rank = -1;
            Node *next = nullptr;
            for(int i = 0; i < N; i++){
				if(current->children[i] != nullptr && I < current->children[i]->I){
                    I = current->children[i]->I;
                    next = current->children[i];
                    road_rank = i;
                }
            }
            if(next == nullptr) return next;
            else{
                int x = current->top[road_rank];
                int y = road_rank;
                board[x][y] = current->player;
                current = next;
            }
        }
		*/
		path.push_back(current);
	}
    return current;
}


int UCT::simulation(Node *node){
    if(node == nullptr) return 1;
    if(node->winner != 0){
        node->Nodewins = 2 * (node->winner == 2 ? 1 : 0);
        return node->Nodewins;
    }
    Node *current = new Node(node->player, N, M, board);
    while(true){
        int choice = 0;
        int setx[M], sety[N];
        for(int i = 0; i < N; i++){
            if(current->top[i] >= 0){
                setx[choice] = current->top[i];
                sety[choice] = i;
                choice++;
            }
        }
        if(choice){
            int x = -1;
            int y = -1;
            for(int i = 0; i < choice && x == -1; i++){
                board[setx[i]][sety[i]] = current->player;
                if(haveWin(setx[i], sety[i])){
                    int result = int(current->player == 2);
                    node->Nodewins = result;
					delete current;
                    return result;
                }
                board[setx[i]][sety[i]] = 0;
            }
            for(int i = 0; i < choice && x == -1; i++){
                board[setx[i]][sety[i]] = (current->player == 1 ? 2 : 1); // nextPlayer
                if(haveWin(setx[i], sety[i])){
                    x = setx[i];
                    y = sety[i];
                }
                board[setx[i]][sety[i]] = 0;
            }
            if(x == -1){
                int t = rand() % choice;
                x = setx[t];
                y = sety[t];
            }
            board[x][y] = current->player;
            if(haveWin(x,y)){
                int result = (current->player == 2 ? 1 : 0);
                node->Nodewins = result;
				delete current;
                return result;
            }
            current->player = (current->player == 1 ? 2 : 1); // nextplayer
            while(x >= 0 && board[x][y] != 0) x--;
            current->top[y] = x;
        }
        else{
            node->Nodewins = 1;
			delete current;
            return 1;
        }
    }
}

void UCT::backPropagation(std::vector<Node*>& path, int result){
    for(int i = 0; i < path.size(); i++){
        path[i]->Nodewins += result;
        path[i]->T += 2;
		if(path[i]->player == 1){ // user
			path[i]->X = double(path[i]->Nodewins) / double(path[i]->T);
			path[i]->I = (root->T <= 1 ? path[i]->X : path[i]->X + sqrt(2 * log(root->T) / double(path[i]->T)));
		}
		else{ // machine
			path[i]->X = double(path[i]->T - path[i]->Nodewins) / double(path[i]->T);
			path[i]->I = (root->T <= 1 ? path[i]->X : path[i]->X + sqrt(2 * log(root->T) / double(path[i]->T)));
		}
    }
}

Point UCT::bestChild(){
	Point result_point(0, 0);
    double X = -100.0;
    for(int i = 0; i < N; i++){
        if(root->children[i] != nullptr && X < root->children[i]->X){
            X = root->children[i]->X;
            result_point.x = root->top[i];
            result_point.y = i;
        }
    }
	return result_point;
}

int UCT::defaultPolicy(Node *node){
    if(node == nullptr) return 1;
    if(node->winner != 0){
        node->Nodewins = 2 * (node->winner == 2 ? 1 : 0);
        return node->Nodewins;
    }
    Node *current = new Node(node->player, N, M, board);
    while(true){
        int choice = 0;
        int setx[M], sety[N];
        for(int i = 0; i < N; i++){
            if(current->top[i] >= 0){
                setx[choice] = current->top[i];
                sety[choice] = i;
                choice++;
            }
        }
        if(choice){
            int x = -1, y = -1;
            for(int i = 0; i < choice && x == -1; i++){
                board[setx[i]][sety[i]] = current->player;
                if(haveWin(setx[i], sety[i])){
                    int result = int(current->player == 2);
                    node->Nodewins = result;
					delete current;
                    return result;
                }
                board[setx[i]][sety[i]] = 0;
            }
            for(int i = 0; i < choice && x == -1; i++){
                board[setx[i]][sety[i]] = (current->player == 1 ? 2 : 1); // nextPlayer
                if(haveWin(setx[i], sety[i])){
                    x = setx[i];
                    y = sety[i];
                }
                board[setx[i]][sety[i]] = 0;
            }
            if(x == -1){
                int t = rand() % choice;
                x = setx[t];
                y = sety[t];
            }
            board[x][y] = current->player;
            if(haveWin(x,y)){
                int result = (current->player == 2 ? 1 : 0);
                node->Nodewins = result;
				delete current;
                return result;
            }
            current->player = (current->player == 1 ? 2 : 1); // nextplayer
            while(x >= 0 && board[x][y] != 0) x--;
            current->top[y] = x;
        }
        else{
            node->Nodewins = 1;
			delete current;
            return 1;
        }
    }
}

void UCT::backup(Node* &pos, std::vector<Node*>& path, int result){
	if(pos){ // 计算该点X、T值
		if(pos->player == 1){ // user
			pos->X = double(pos->Nodewins) / double(pos->T);
			pos->I = (root->T <= 1 ? pos->X : pos->X + sqrt(2 * log(root->T) / double(pos->T)));
		}
		else{ // machine
			pos->X = double(pos->T - pos->Nodewins) / double(pos->T);
			pos->I = (root->T <= 1 ? pos->X : pos->X + sqrt(2 * log(root->T) / double(pos->T)));
		}
	} 
    for(int i = 0; i < path.size(); i++){
        path[i]->Nodewins += result;
        path[i]->T += 2;
		if(path[i]->player == 1){ // user
			path[i]->X = double(path[i]->Nodewins) / double(path[i]->T);
			path[i]->I = (root->T <= 1 ? path[i]->X : path[i]->X + sqrt(2 * log(root->T) / double(path[i]->T)));
		}
		else{ // machine
			path[i]->X = double(path[i]->T - path[i]->Nodewins) / double(path[i]->T);
			path[i]->I = (root->T <= 1 ? path[i]->X : path[i]->X + sqrt(2 * log(root->T) / double(path[i]->T)));
		}
    }
}

Point UCT::getPoint(){
	while((double)(clock() - startTime) / CLOCKS_PER_SEC < TIME){ // 在时间限制内
		for(int i = 0; i < M; i++) for(int j = 0; j < N; j++) board[i][j] = baseBoard[i][j]; // 复原棋盘
		std::vector<Node*> path;
		Node *pos = treePolicy(path); // 选择一个点并记录路径
		
		int result = 0;
		if((pos && pos->winner == 0) || !pos) result = defaultPolicy(pos); // 如果选到非终止节点，则模拟
		else result = 2 * (pos->winner == 2 ? 1 : 0); // 如果选到终止节点，直接更新
		
		backup(pos, path, result); // 回溯更新该路径
	}

	/*
	struct timeval currentTime; // 用于计时
	double timeInterval = 0.0; // 单位为微妙
    while(timeInterval < TIME){ // 时间不超过TIME的情况下进行while循环
		for(int i = 0; i < M; i++){
            for(int j = 0; j < N; j++) board[i][j] = baseBoard[i][j]; // 复原棋盘
		}
		std::vector<Node*> path;
		Node *pos = selection(path); // 选择一个点并记录路径
		int result = 0;
		if((pos && pos->winner == 0) || !pos) result = simulation(pos); // 如果选到非终止节点，则模拟
		else result = 2 * (pos->winner == 2 ? 1 : 0); // 如果选到终止节点，直接更新
		if(pos){ // 计算该点Xj、Tj值
			if(pos->player == 1){ // user
				pos->X = double(pos->Nodewins) / double(pos->T);
				pos->I = (root->T <= 1 ? pos->X : pos->X + sqrt(2 * log(root->T) / double(pos->T)));
			}
			else{ // machine
				pos->X = double(pos->T - pos->Nodewins) / double(pos->T);
				pos->I = (root->T <= 1 ? pos->X : pos->X + sqrt(2 * log(root->T) / double(pos->T)));
			}
		} 
		backPropagation(path, result); // 回溯更新该路径
		gettimeofday(&currentTime, NULL); // 更新当前时间
		timeInterval = (currentTime.tv_sec - startTime.tv_sec) * 1000000 + (currentTime.tv_usec - startTime.tv_usec); // 单位为微秒
	}
	*/

	Point best_child = bestChild();
	return best_child;
	/*
	Point result_point(0, 0);
    double X = -100.0;
    for(int i = 0; i < N; i++){
        if(root->children[i] != nullptr && X < root->children[i]->X){
            X = root->children[i]->X;
            result_point.x = root->top[i];
            result_point.y = i;
        }
    }
	return result_point;
	*/
}