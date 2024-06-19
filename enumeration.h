#pragma once
#include "graph.h"
#include "BuildTable.h"

int isolated_localtion; 
bool is_isolated[260];
int reverse_index[260];
ui *embedding;
int parent[260];
vector <int> act_isolates[260];

bool same2next[256];
int isolated_number = 0;
BSRGraph*** qfliter_bsr_graph_;

void allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                              ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui **&valid_candidate_idx) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num) {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        valid_candidate_idx[i] = new ui[max_candidates_num];
        embedding[i] = -1;
    }
}

void releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
    ui **valid_candidate_idx) {
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    for (ui i = 0; i < query_vertices_num; ++i) {
        delete[] valid_candidate_idx[i];
    }

    delete[] valid_candidate_idx;
}

std::unordered_set<VertexID> vis_right;
std::unordered_map<VertexID, VertexID> reverse_embedding;
vector <int> fail_vec; 
bool find_augmentation(const int &lim, int cur_u, ui **candidates, ui *idx_count, ui **valid_candidate_idx) {
    int u_rank = reverse_index[cur_u];
    fail_vec.push_back(cur_u);
    if (u_rank <= lim ) return false;
    assert(reverse_index[cur_u] != lim);

    for (int i = 0; i < idx_count[u_rank]; ++i) {
        int idx = valid_candidate_idx[u_rank][i];
        int v = candidates[cur_u][idx];
        if (vis_right.find(v) != vis_right.end())
            continue;
        vis_right.insert(v);
        if (reverse_embedding.find(v) == reverse_embedding.end() || find_augmentation(lim, reverse_embedding[v], candidates, idx_count, valid_candidate_idx)) {
            reverse_embedding[v] = cur_u;
            embedding[cur_u] = v;
            return true;
        }
    }
    return false;
}

long long tmp_edge = 0;
long long ug_edges = 0;
long long ug_num = 0;
std::bitset<260> vg_ancestors[260][260*2];
BSRSet bsrs[260][260*2];
int updated_times[260];
vector<int> update_vertex[260];

bool in_set(int ele, ui *candidates, int cnt, ui *valid_candidate_index) {
    for (int j = 0; j < cnt; ++j) {
        int idx = valid_candidate_index[j];
        int v = candidates[idx];

        if (v == ele) 
            return true;
        if (v > ele) 
            return false;
    }
    return false;
}

int generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index, ui **an, ui *an_cnt, ui *order, ui *candidates_count, ui **candidates, const Graph* query_graph) {
    VertexID u = order[depth];
    ui current_index_id = idx_embedding[u];
    update_vertex[depth].clear();
    for (ui i = 0; i < an_cnt[depth]; ++i) {
        VertexID current_an = an[depth][i];
        int an_idx = reverse_index[current_an];
        int t = updated_times[current_an];

        BSRGraph &cur_bsr_graph = *qfliter_bsr_graph_[u][current_an];
        BSRSet &cur_bsr_set = cur_bsr_graph.bsrs[current_index_id];
        if (bsrs[current_an][t].base_ == nullptr) {
            bsrs[current_an][t].base_ = new int[candidates_count[current_an]];
            bsrs[current_an][t].states_ = new int[candidates_count[current_an]];
        }
        if (t > 0) {
            bsrs[current_an][t].size_ = intersect_qfilter_bsr_b4_v2(cur_bsr_set.base_, cur_bsr_set.states_,
                cur_bsr_set.size_, bsrs[current_an][t-1].base_, bsrs[current_an][t-1].states_, bsrs[current_an][t-1].size_, bsrs[current_an][t].base_, bsrs[current_an][t].states_);

            if (bsrs[current_an][t-1] != bsrs[current_an][t]) {
                vg_ancestors[current_an][t+1] = vg_ancestors[current_an][t] | vg_ancestors[u][updated_times[u]];
                update_vertex[depth].push_back(current_an);
            } else {
                updated_times[current_an]--;
            }
        } else {
            bsrs[current_an][t].size_ = cur_bsr_set.size_;
            memcpy(bsrs[current_an][t].base_, cur_bsr_set.base_, sizeof(int) * cur_bsr_set.size_);
            memcpy(bsrs[current_an][t].states_, cur_bsr_set.states_, sizeof(int) * cur_bsr_set.size_);
            vg_ancestors[current_an][t+1] = vg_ancestors[current_an][t] | vg_ancestors[u][updated_times[u]];
            update_vertex[depth].push_back(current_an);
        }
        updated_times[current_an]++;
        if (bsrs[current_an][t].size_ == 0)
            return current_an;

        idx_count[an_idx] = offline_bsr_trans_uint(bsrs[current_an][t].base_, bsrs[current_an][t].states_, bsrs[current_an][t].size_, (int *) valid_candidate_index[an_idx]);
        if (embedding[current_an] != -1) {
            if (in_set(embedding[current_an], candidates[current_an], idx_count[an_idx], valid_candidate_index[an_idx])) {
                continue;
            } else {
                reverse_embedding.erase(embedding[current_an]);
                embedding[current_an] = -1;
            }
        }
        vis_right.clear();
        bool res = find_augmentation(depth, current_an, candidates, idx_count, valid_candidate_index);
        if (res) {
            fail_vec.clear();
        } else {
            for (int d = 0; d < i ; ++d) {
                int current_an = an[depth][d];
                if (updated_times[current_an] > 1) {
                    reverse_embedding.erase(embedding[current_an]);
                    embedding[current_an] = -1;
                }
            }
            for (auto &un: fail_vec){
                vg_ancestors[current_an][updated_times[current_an]] |= vg_ancestors[un][updated_times[un]];
            }
            fail_vec.clear();
            return current_an;
        }
    }
    return -1;
}

int perm(int level, ui *order, int max_depth, int output_limit_num, ui **candidates, ui *idx_count, ui **valid_candidate_idx, size_t &c) {
    if (level == max_depth) {
        return 1;
    }

    int u = order[level];
    assert(is_isolated[u]); 
    int ans = 0;
    int delta = -1;

    for (int i = 0; i < idx_count[level]; ++i) {
        int idx = valid_candidate_idx[level][i]; 
        int v = candidates[u][idx]; 
        assert(reverse_embedding.find(embedding[u]) != reverse_embedding.end());
        reverse_embedding.erase(embedding[u]);
        vis_right.clear();
        vis_right.insert(v);
        if (reverse_embedding.find(v) != reverse_embedding.end()) {
            assert(reverse_embedding[v] != u); 
        }
        if (reverse_embedding.find(v) == reverse_embedding.end() || find_augmentation(level, reverse_embedding[v], candidates, idx_count, valid_candidate_idx)) {
            reverse_embedding[v] = u;
            embedding[u] = v;
            if (!same2next[level]) {
                if (delta == -1) 
                    delta = perm(level + 1, order, max_depth, output_limit_num, candidates, idx_count, valid_candidate_idx, c);
                assert(delta != -1);
                ans += delta;
            } else
                ans += perm(level + 1, order, max_depth, output_limit_num, candidates, idx_count, valid_candidate_idx, c);
            if (c + ans >= output_limit_num) 
                return ans;
        } else {
            reverse_embedding[embedding[u]] = u;
        }
    }
    return ans;
}

int fa[6000000];
bool visited[260];
ui new_order[260];
ui real_level[260];
vector<int> visited_data_vertices;
int find(int x) {
    if (fa[x] == x) 
        return x;
    return fa[x] = find(fa[x]);
}
void counting(const Graph *query_graph, ui *order, int output_limit_num, ui **candidates, ui *candidates_count, ui *idx_count, ui **valid_candidate_idx, size_t &c) {
    int query_vertex_number = query_graph->getVerticesCount();
    assert(reverse_embedding.size() == query_vertex_number);
    int ptr = isolated_localtion;

    // construct new_order by connected component
    for (int i = isolated_localtion; i < query_vertex_number; ++i) {
        int u = order[i];
        fa[u] = u;
        for (int j = 0; j < idx_count[i]; ++j) {
            int idx = valid_candidate_idx[i][j];
            int v = candidates[order[i]][idx];
            if (fa[v+query_vertex_number] == -1){
                visited_data_vertices.push_back(v);
                fa[v+query_vertex_number] = find(u);
            } else if (find(v+query_vertex_number) != find(u)) {
                fa[find(v+query_vertex_number)] = find(u);
            }
        }
    }
    for (int i = isolated_localtion; i < query_vertex_number; ++i) if (!visited[order[i]]) {
        visited[order[i]] = true;
        same2next[ptr] = true;
        real_level[ptr] = i;
        new_order[ptr++] = order[i];
        for (int j = i + 1; j < query_vertex_number; ++j) {
            if (find(order[i]) == find(order[j])) {
                visited[order[j]] = true;
                same2next[ptr] = true;
                real_level[ptr] = j;
                new_order[ptr++] = order[j];
            }
        }
        same2next[ptr - 1] = false;
    }
    assert(ptr == query_vertex_number);
    ui tmp_idx_count[260];
    ui *tmp_valid_candidate_idx[260];
    // change idx_count and valid_candidate_idx by real_level
    for (int i = isolated_localtion; i < query_vertex_number; ++i) {
        tmp_idx_count[i] = idx_count[real_level[i]];
        tmp_valid_candidate_idx[i] = valid_candidate_idx[real_level[i]];
    }
    for (int i = isolated_localtion; i < query_vertex_number; ++i) {
        idx_count[i] = tmp_idx_count[i];
        valid_candidate_idx[i] = tmp_valid_candidate_idx[i];
        order[i] = new_order[i];
        reverse_index[order[i]] = i;
    }

    c += perm(isolated_localtion, order, query_vertex_number, output_limit_num, candidates, idx_count, valid_candidate_idx, c);
    if (c > output_limit_num) {
        c = output_limit_num;
    }

    for (auto &v: visited_data_vertices) 
        fa[v + query_vertex_number] = -1;
    visited_data_vertices.clear();
    for (int i = isolated_localtion; i < query_vertex_number; ++i) 
        visited[order[i]] = false;
}

void generateAN(const Graph *query_graph, ui *order, ui **&an, ui *&an_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    an_count = new ui[query_vertices_num];
    std::fill(an_count, an_count + query_vertices_num, 0);
    an = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        an[i] = new ui[query_vertices_num];
        parent[i] = -1;
    }
    std::vector<bool> visited_vertices(query_vertices_num, false);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        ui nbrs_cnt;
        assert(reverse_index[vertex] == i);
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];
            if (!visited_vertices[nbr]) {
                an[i][an_count[i]++] = nbr;
                if (parent[nbr] == -1) 
                    parent[nbr] = vertex;
            }
        }
        visited_vertices[vertex] = true;
    }
    assert(parent[order[0]] == -1);
    for (ui i = 1; i < query_vertices_num; ++i) {
        assert(parent[order[i]] != -1);
        act_isolates[parent[order[i]]].push_back(order[i]);
    }
}

void identifying(const Graph *query_graph, ui *order, ui* candidates_count) {
    bool used[256];
    int N = query_graph->getVerticesCount();
    for (int i = 0; i < N; ++i) {
        is_isolated[i] = false;
        used[i] = false;
    }
    ui edge_count;
    for (int loop = 0; loop < N; ++loop) {
        int best = order[loop];
        assert(!used[best]);
        const ui* neighbors = query_graph -> getVertexNeighbors(best, edge_count);
        is_isolated[best] = true;
        for (int p = 0; p < edge_count; ++p) {
            int un = neighbors[p];
            if (!used[un]) {
                is_isolated[best] = false;
                break;
            }
        }
        if (!is_isolated[best])
            isolated_localtion++;
        used[best] = true;
    }

    for (int i = 0; i < N; ++i) {
        int u = order[i];
        reverse_index[u] = i;
    }
}

size_t explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, 
    ui *candidates_count, size_t output_limit_num, size_t &call_count, ui* matching_order) {
    ui *order = matching_order;
    identifying(query_graph, order, candidates_count);

    // Generate an.
    ui **an;
    ui *an_count;
    generateAN(query_graph, order, an, an_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *idx_embedding;
    ui **valid_candidate_idx;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                valid_candidate_idx);
    memset(fa, -1, sizeof(int) * 6000000);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<260>> vec_failing_set(max_depth);
    reverse_embedding.reserve(260 * 2);
    for (int i = 0; i < max_depth; ++i) {
        vg_ancestors[i][0].reset();
        vg_ancestors[i][0].set(i);
        embedding[i] = -1;
        if (i != 0) idx_count[i] = 0;
    }
    double avg_failing_set = 0;
    long long call_failing_set = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            call_count += 1;

            if (cur_depth == isolated_localtion ) {
                counting(query_graph, order, output_limit_num, candidates, candidates_count, idx_count, valid_candidate_idx, embedding_cnt);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
                break;
            }

            if (embedding[u] != -1) {
                reverse_embedding.erase(embedding[u]);
                embedding[u] = -1;
            }

            if (reverse_embedding.find(v) != reverse_embedding.end()) {
                vis_right.clear();
                vis_right.insert(v);
                if (find_augmentation(cur_depth, reverse_embedding[v], candidates, idx_count, valid_candidate_idx)) {
                    fail_vec.clear();
                } else {
                    idx[cur_depth] += 1;
                    vec_failing_set[cur_depth] = vg_ancestors[u][updated_times[u]];
                    for (int i=0;i< fail_vec.size();++i) {
                        int un = fail_vec[i];
                        vec_failing_set[cur_depth] |= vg_ancestors[un][updated_times[un]];
                    }
                    fail_vec.clear();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                    continue;
                }
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            reverse_embedding[v] = u;
            idx[cur_depth] += 1;

            int flag = generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, an, an_count, order, candidates_count, candidates, query_graph);
            cur_depth += 1;
            idx[cur_depth] = 0;

            if (flag >= 0) {
                vec_failing_set[cur_depth - 1] = vg_ancestors[flag][updated_times[flag]];
                break;
            } else {
                vec_failing_set[cur_depth - 1].reset();
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
            for (auto an: update_vertex[cur_depth]) {
                --updated_times[an];
                int an_idx = reverse_index[an];
                int t = updated_times[an] - 1;
                if (t >= 0) 
                    idx_count[an_idx] = offline_bsr_trans_uint(bsrs[an][t].base_, bsrs[an][t].states_, bsrs[an][t].size_, (int *) valid_candidate_idx[an_idx]);
            }

            for (int i = 0; i < act_isolates[u].size(); ++i) {
                int an = act_isolates[u][i];
                reverse_embedding.erase(embedding[an]);
                embedding[an] = -1;
                idx_count[reverse_index[an]] = 0;
            }

            if (cur_depth != 0) { 
                int num = 0;
                for (int i = 0;i < cur_depth; ++i) {
                    if (vec_failing_set[cur_depth].test(i)) {
                        num += 1;
                    }
                }
                avg_failing_set += num * 1. / (cur_depth);
                ++call_failing_set;
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
        }
    }
    // Release the buffer.
    EXIT:

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    printf("Enumeration Time(second): %f\n", duration.count());

    return embedding_cnt;
}
