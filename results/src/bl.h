//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include "parallel.h"
#include "init.h"
class bl{
public:
	void bl_init(block_t *bb, std::vector<int> &ijk_sg);
	void bl_clear_all_ghosts(bool old);
	void bl_remove_from_list(block_t **bb, block_t **bm, bool first);
	void bl_free(block_t *bb);
	void bl_activate_new(block_t *bb, std::vector<int> &ijk_sg);
};
