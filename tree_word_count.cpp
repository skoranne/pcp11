////////////////////////////////////////////////////////////////////////////////
// File     : tree_word_count.cpp
// Author   : Sandeep Koranne (C) 2021. All rights reserved.
// Purpose  : Count occurences of words using Tree
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <cassert>
#include <vector>

namespace TreeCount {

  struct Node {
    std::string word;
    int count;
    Node *left=nullptr, *right=nullptr;
    Node( const std::string& inword ): word{ inword }, count{1} {}
  };
  using V = std::vector<std::string>;
  Node* AddTree( Node* p, const std::string& word ) {
    if( p == nullptr ) {
      p = new Node( word );
      return p;
    }
    if( word == p->word ) { p->count++; return p; }
    if( word < p->word ) 
      p->left = AddTree( p->left, word );
    else
      p->right = AddTree( p->right, word );
    return p;
  }
  // Inorder traversal
  void PrintTree( Node* p ) {
    if( p->left ) PrintTree( p->left );
    std::cout << p->word << "\t" << p->count << std::endl;
    if( p->right ) PrintTree( p->right );
  }
    
  void CountWords( const V& words ) {
    if( words.empty() ) return;
    Node * root = nullptr;
    for( auto x : words ) root = AddTree( root, x );
    assert( root && "Root should have been built" );
    PrintTree( root );
  }
}
using namespace TreeCount;
int main() {
  V word_vec;
  while( std::cin ) {
    std::string word;
    std::cin >> word;
    if( word != "" ) word_vec.push_back( word );
  }
  //for( auto x : word_vec ) std::cout << x << " ";
  CountWords( word_vec );
  return 0;
}

