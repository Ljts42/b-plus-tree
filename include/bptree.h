#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace tree_details {

template <class T>
constexpr bool IsConst = false;
template <template <bool> class T>
constexpr bool IsConst<T<true>> = true;

} // namespace tree_details

template <class Key, class Value, std::size_t BlockSize = 4096, class Less = std::less<Key>>
class BPTree
{
    struct Node;
    struct Leaf;
    struct Inner;

    template <bool is_const>
    class Iterator;

public:
    using key_type = Key;
    using mapped_type = Value;
    using value_type = std::pair<Key, Value>; // NB: a digression from std::map
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using size_type = std::size_t;

    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;

    BPTree();
    BPTree(std::initializer_list<std::pair<Key, Value>>);
    BPTree(const BPTree &);
    BPTree(BPTree &&);

    BPTree & operator=(const BPTree &);
    BPTree & operator=(BPTree &&);

    iterator begin();
    const_iterator cbegin() const;
    const_iterator begin() const;
    iterator end();
    const_iterator cend() const;
    const_iterator end() const;

    bool empty() const;
    size_type size() const;
    void clear();

    size_type count(const Key &) const;
    bool contains(const Key &) const;
    std::pair<iterator, iterator> equal_range(const Key &);
    std::pair<const_iterator, const_iterator> equal_range(const Key &) const;
    iterator lower_bound(const Key &);
    const_iterator lower_bound(const Key &) const;
    iterator upper_bound(const Key &);
    const_iterator upper_bound(const Key &) const;
    iterator find(const Key & key);
    const_iterator find(const Key & key) const;

    // 'at' method throws std::out_of_range if there is no such key
    Value & at(const Key &);
    const Value & at(const Key &) const;

    // '[]' operator inserts a new element if there is no such key
    Value & operator[](const Key &);

    std::pair<iterator, bool> insert(const Key &, const Value &); // NB: a digression from std::map
    std::pair<iterator, bool> insert(const Key &, Value &&);      // NB: a digression from std::map
    template <class ForwardIt>
    void insert(ForwardIt, ForwardIt);
    void insert(std::initializer_list<value_type>);
    iterator erase(const_iterator);
    iterator erase(const_iterator, const_iterator);
    size_type erase(const Key &);

    ~BPTree();

private:
    static constexpr size_type m_order = std::max<std::size_t>(2, (BlockSize - sizeof(std::size_t) - sizeof(void *)) / (sizeof(Key) + sizeof(void *)));
    static constexpr Less less{};
    Node * m_root = nullptr;
    size_type m_tree_size = 0;

    static Leaf * find_min_leaf(Node *);
    static Leaf * find_max_leaf(Node *);

    static Node * copy(Node *);
    void redirect();

    static void clear(Node *);

    std::pair<Leaf *, size_type> search(const Key &) const;

    template <class K, class V>
    iterator insert(Leaf *, size_type, K &&, V &&);

    void split(Leaf *);
    void split(Inner *);

    static Leaf * find_left(Leaf *);
    static Inner * find_left(Inner *);
    static Leaf * find_right(Leaf *);
    static Inner * find_right(Inner *);

    bool erase(Leaf *, size_type, const Key &);
    static void erase(Inner *, size_type, const Key &);

    static void merge(Leaf *, Leaf *, const Key &);
    static void merge(Inner *, Inner *, const Key &);

    static void update(Node *, const Key &, const Key &);
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Node
{
    Inner * m_parent = nullptr;
    size_type m_size = 0;

    Node() = default;
    Node(const Node & node)
        : m_size(node.m_size)
    {
    }

    virtual size_type size() const
    {
        return m_size;
    }

    virtual const Key & get_key(size_type) const = 0;
    virtual size_type get_ind(const Key &) const = 0;
    virtual Node * get_child(size_type) const = 0;

    virtual void shift(size_type) = 0;
    virtual void erase(size_type, size_type) = 0;

    virtual ~Node() = default;
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Leaf : Node
{
    Leaf * m_right = nullptr;
    std::array<value_type, m_order> m_data{};

    Leaf() = default;
    Leaf(const Leaf & node)
        : Node(node)
        , m_data(node.m_data)
    {
    }

    const Key & get_key(size_type ind) const override
    {
        return m_data[ind].first;
    }

    Node * get_child(size_type) const override
    {
        return nullptr;
    }

    size_type get_ind(const Key & key) const override
    {
        auto compare = [](const auto & lhs, const auto & rhs) {
            return less(lhs.first, rhs);
        };
        return std::lower_bound(m_data.begin(), m_data.begin() + this->m_size, key, compare) - m_data.begin();
    }

    void shift(size_type begin) override
    {
        for (size_type i = this->m_size; i > begin; --i) {
            std::swap(m_data[i - 1], m_data[i]);
        }
        ++this->m_size;
    }

    void erase(size_type begin, size_type end) override
    {
        for (size_type i = 0; end + i < this->m_size; ++i) {
            std::swap(m_data[begin + i], m_data[end + i]);
        }
        this->m_size -= end - begin;
    }
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Inner : Node
{
    std::array<key_type, m_order> m_keys{};
    std::array<Node *, m_order + 1> m_children{};

    Inner() = default;
    Inner(const Inner & node)
        : Node(node)
        , m_keys(node.m_keys)
    {
    }

    const Key & get_key(size_type ind) const override
    {
        return m_keys[ind];
    }

    Node * get_child(size_type ind) const override
    {
        return (ind <= this->m_size) ? m_children[ind] : nullptr;
    }

    size_type get_ind(const Key & key) const override
    {
        return std::upper_bound(m_keys.begin(), m_keys.begin() + this->m_size, key, less) - m_keys.begin();
    }

    void shift(size_type begin) override
    {
        for (size_type i = this->m_size; i > begin; --i) {
            std::swap(m_keys[i - 1], m_keys[i]);
            std::swap(m_children[i], m_children[i + 1]);
        }
        ++this->m_size;
    }

    void erase(size_type begin, size_type end) override
    {
        for (size_type i = 0; i + end < this->m_size; ++i) {
            std::swap(m_keys[begin + i], m_keys[end + i]);
            std::swap(m_children[begin + i + 1], m_children[end + i + 1]);
        }
        for (size_type i = begin; i < end; ++i) {
            m_children[this->m_size] = nullptr;
            --this->m_size;
        }
    }
};

template <class Key, class Value, std::size_t BlockSize, class Less>
template <bool is_const>
class BPTree<Key, Value, BlockSize, Less>::Iterator
{
public:
    using difference_type = std::ptrdiff_t;
    using value_type = std::conditional_t<is_const, const std::pair<Key, Value>, std::pair<Key, Value>>;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::forward_iterator_tag;

    Iterator() = default;

    template <class R = Iterator, std::enable_if_t<tree_details::IsConst<R>, int> = 0>
    Iterator(const Iterator<false> & other)
        : m_current(other.m_current)
        , m_pos(other.m_pos)
    {
    }

    friend bool operator==(const Iterator & lhs, const Iterator & rhs)
    {
        return lhs.m_current == rhs.m_current && lhs.m_pos == rhs.m_pos;
    }

    friend bool operator!=(const Iterator & lhs, const Iterator & rhs)
    {
        return !(lhs == rhs);
    }

    reference operator*() const
    {
        return m_current->m_data[m_pos];
    }

    pointer operator->() const
    {
        return &m_current->m_data[m_pos];
    }

    Iterator & operator++()
    {
        ++m_pos;
        if (m_pos >= m_current->size() && m_current->m_right) {
            m_pos = 0;
            m_current = m_current->m_right;
        }
        return *this;
    }

    Iterator operator++(int)
    {
        auto tmp = *this;
        operator++();
        return tmp;
    }

private:
    friend class BPTree;

    Iterator(Leaf * node)
        : m_current(node)
    {
    }

    Iterator(Leaf * node, const size_type & pos)
        : m_current(node)
        , m_pos(pos)
    {
    }

    Leaf * m_current = nullptr;
    size_type m_pos = 0;
};

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree()
    : m_root(new Leaf)
{
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(std::initializer_list<std::pair<Key, Value>> list)
    : m_root(new Leaf)
{
    insert(list);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(const BPTree & tree)
    : m_root(copy(tree.m_root))
    , m_tree_size(tree.m_tree_size)
{
    redirect();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::operator=(const BPTree & tree)
        -> BPTree &
{
    clear(m_root);
    m_root = copy(tree.m_root);
    m_tree_size = tree.m_tree_size;
    redirect();
    return *this;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::copy(Node * node)
        -> Node *
{
    if (!node->get_child(0)) {
        return new Leaf(*dynamic_cast<Leaf *>(node));
    }
    else {
        Inner * inner = new Inner(*dynamic_cast<Inner *>(node));
        for (size_type i = 0; i <= node->size(); ++i) {
            Node * child = copy(node->get_child(i));
            inner->m_children[i] = child;
            inner->m_children[i]->m_parent = inner;
        }
        return inner;
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::redirect()
{
    if (m_root && m_root->get_child(0)) {
        Leaf * current_leaf = find_min_leaf(m_root);
        Inner * curent_inner = current_leaf->m_parent;

        for (size_type ind = 1; ind <= curent_inner->size(); ++ind) {
            Leaf * right_leaf = dynamic_cast<Leaf *>(curent_inner->get_child(ind));
            current_leaf->m_right = right_leaf;
            current_leaf = right_leaf;
        }

        Inner * right_inner = find_right(curent_inner);
        while (right_inner) {
            curent_inner = right_inner;
            current_leaf->m_right = dynamic_cast<Leaf *>(curent_inner->get_child(0));
            current_leaf = current_leaf->m_right;

            for (size_type ind = 1; ind <= curent_inner->size(); ++ind) {
                Leaf * right_leaf = dynamic_cast<Leaf *>(curent_inner->get_child(ind));
                current_leaf->m_right = right_leaf;
                current_leaf = right_leaf;
            }
            right_inner = find_right(curent_inner);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(BPTree && tree)
    : m_root(new Leaf)
{
    std::swap(m_root, tree.m_root);
    std::swap(m_tree_size, tree.m_tree_size);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::operator=(BPTree && tree)
        -> BPTree &
{
    std::swap(m_root, tree.m_root);
    std::swap(m_tree_size, tree.m_tree_size);
    return *this;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_min_leaf(Node * node)
        -> Leaf *
{
    while (node->get_child(0)) {
        node = node->get_child(0);
    }
    return dynamic_cast<Leaf *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_max_leaf(Node * node)
        -> Leaf *
{
    while (node->get_child(node->size())) {
        node = node->get_child(node->size());
    }
    return dynamic_cast<Leaf *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::begin()
        -> iterator
{
    return iterator(find_min_leaf(m_root));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::cbegin() const
        -> const_iterator
{
    return const_iterator(find_min_leaf(m_root));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::begin() const
        -> const_iterator
{
    return cbegin();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::end()
        -> iterator
{
    Leaf * leaf = find_max_leaf(m_root);
    return iterator(leaf, leaf->size());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::cend() const
        -> const_iterator
{
    Leaf * leaf = find_max_leaf(m_root);
    return const_iterator(leaf, leaf->size());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::end() const
        -> const_iterator
{
    return cend();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::empty() const
{
    return size() == 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::size() const
        -> size_type
{
    return m_tree_size;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::count(const Key & key) const
        -> size_type
{
    return contains(key) ? 1 : 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::contains(const Key & key) const
{
    const_iterator lower = lower_bound(key);
    return lower != end() && !less(key, lower->first);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::equal_range(const Key & key)
        -> std::pair<iterator, iterator>
{
    return {lower_bound(key), upper_bound(key)};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::equal_range(const Key & key) const
        -> std::pair<const_iterator, const_iterator>
{
    return {lower_bound(key), upper_bound(key)};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::search(const Key & key) const
        -> std::pair<Leaf *, size_type>
{
    Node * node = m_root;
    while (node->get_child(0)) {
        size_type ind = node->get_ind(key);
        node = node->get_child(ind);
    }
    Leaf * leaf = dynamic_cast<Leaf *>(node);
    size_type ind = leaf->get_ind(key);
    return {leaf, ind};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::lower_bound(const Key & key)
        -> iterator
{
    auto [leaf, ind] = search(key);
    if (ind >= leaf->size() && leaf->m_right) {
        return iterator(leaf->m_right);
    }
    return iterator(leaf, ind);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::lower_bound(const Key & key) const
        -> const_iterator
{
    auto [leaf, ind] = search(key);
    if (ind >= leaf->size() && leaf->m_right) {
        return const_iterator(leaf->m_right);
    }
    return const_iterator(leaf, ind);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::upper_bound(const Key & key)
        -> iterator
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        ++lower;
    }
    return lower;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::upper_bound(const Key & key) const
        -> const_iterator
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        ++lower;
    }
    return lower;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find(const Key & key)
        -> iterator
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        return lower;
    }
    return end();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find(const Key & key) const
        -> const_iterator
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        return lower;
    }
    return cend();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline Value & BPTree<Key, Value, BlockSize, Less>::at(const Key & key)
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        return lower->second;
    }
    throw std::out_of_range("");
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline const Value & BPTree<Key, Value, BlockSize, Less>::at(const Key & key) const
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        return lower->second;
    }
    throw std::out_of_range("");
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline Value & BPTree<Key, Value, BlockSize, Less>::operator[](const Key & key)
{
    auto [leaf, ind] = search(key);
    iterator it = iterator(leaf, ind);
    if (ind >= leaf->size() || less(key, leaf->get_key(ind))) {
        it = insert(leaf, ind, key, Value());
        ++m_tree_size;
    }
    return it->second;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(const Key & key, const Value & value)
        -> std::pair<iterator, bool>
{
    auto [leaf, ind] = search(key);
    if (ind >= leaf->size() || less(key, leaf->get_key(ind))) {
        iterator it = insert(leaf, ind, key, value);
        ++m_tree_size;
        return {it, true};
    }
    return {iterator(leaf, ind), false};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(const Key & key, Value && value)
        -> std::pair<iterator, bool>
{
    auto [leaf, ind] = search(key);
    if (ind >= leaf->size() || less(key, leaf->get_key(ind))) {
        iterator it = insert(leaf, ind, key, std::forward<Value>(value));
        ++m_tree_size;
        return {it, true};
    }
    return {iterator(leaf, ind), false};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
template <class ForwardIt>
inline void BPTree<Key, Value, BlockSize, Less>::insert(ForwardIt begin, ForwardIt end)
{
    for (auto i = begin; i != end; ++i) {
        insert(i->first, i->second);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::insert(std::initializer_list<value_type> list)
{
    for (const auto & i : list) {
        insert(i.first, i.second);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
template <class K, class V>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(Leaf * leaf, size_type ind, K && key, V && value)
        -> iterator
{
    leaf->shift(ind);
    leaf->m_data[ind] = {std::forward<K>(key), std::forward<V>(value)};

    if (leaf->size() == m_order) {
        split(leaf);
        if (ind >= m_order / 2) {
            ind -= m_order / 2;
            leaf = leaf->m_right;
        }
    }
    return iterator(leaf, ind);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::split(Leaf * left)
{
    Leaf * right = new Leaf;
    right->m_right = left->m_right;
    left->m_right = right;

    for (size_type i = m_order / 2; i < left->size(); ++i) {
        std::swap(left->m_data[i], right->m_data[right->size()]);
        ++right->m_size;
    }
    left->erase(m_order / 2, left->size());

    if (!left->m_parent) {
        Inner * parent = new Inner;
        parent->m_keys[0] = right->get_key(0);
        parent->m_children[0] = left;
        parent->m_children[0]->m_parent = parent;
        parent->m_children[1] = right;
        parent->m_children[1]->m_parent = parent;
        ++parent->m_size;
        m_root = parent;
    }
    else {
        Inner * parent = left->m_parent;
        size_type ind = std::lower_bound(parent->m_keys.begin(), parent->m_keys.begin() + parent->size(), right->get_key(0), less) - parent->m_keys.begin();
        parent->shift(ind);
        parent->m_keys[ind] = right->get_key(0);
        parent->m_children[ind + 1] = right;
        parent->m_children[ind + 1]->m_parent = parent;
        if (parent->size() == m_order) {
            split(parent);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::split(Inner * left)
{
    Inner * right = new Inner;
    for (size_type i = m_order / 2 + 1; i < left->size(); ++i) {
        std::swap(right->m_keys[right->size()], left->m_keys[i]);
        std::swap(right->m_children[right->size()], left->m_children[i]);
        right->m_children[right->size()]->m_parent = right;
        ++right->m_size;
    }
    std::swap(right->m_children[right->size()], left->m_children[left->size()]);
    right->m_children[right->size()]->m_parent = right;

    left->erase(m_order / 2, left->size());

    if (!left->m_parent) {
        Inner * parent = new Inner;
        parent->m_keys[0] = find_min_leaf(right)->get_key(0);
        parent->m_children[0] = left;
        parent->m_children[0]->m_parent = parent;
        parent->m_children[1] = right;
        parent->m_children[1]->m_parent = parent;
        ++parent->m_size;
        m_root = parent;
    }
    else {
        Inner * parent = left->m_parent;
        Leaf * right_min = find_min_leaf(right);
        size_type ind = std::lower_bound(parent->m_keys.begin(), parent->m_keys.begin() + parent->size(), right_min->get_key(0), less) - parent->m_keys.begin();
        parent->shift(ind);
        parent->m_keys[ind] = right_min->get_key(0);
        parent->m_children[ind + 1] = right;
        parent->m_children[ind + 1]->m_parent = parent;
        if (parent->size() == m_order) {
            split(parent);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Leaf * leaf)
        -> Leaf *
{
    if (!leaf || !leaf->m_parent) {
        return nullptr;
    }
    Inner * parent = leaf->m_parent;
    size_type ind = parent->get_ind(leaf->get_key(0));
    if (ind != 0) {
        return dynamic_cast<Leaf *>(parent->get_child(ind - 1));
    }
    Inner * left = find_left(parent);
    if (!left) {
        return nullptr;
    }
    return dynamic_cast<Leaf *>(left->get_child(left->size()));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Leaf * node)
        -> Leaf *
{
    return node->m_right;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Inner * inner)
        -> Inner *
{
    if (!inner || !inner->m_parent) {
        return nullptr;
    }
    const Key & key = inner->get_key(0);
    Inner * parent = inner->m_parent;
    size_type height = 0;
    size_type ind = 0;
    while (parent && ind == 0) {
        inner = parent;
        ++height;
        parent = inner->m_parent;
        ind = inner->get_ind(key);
    }
    if (!parent && ind == 0) {
        return nullptr;
    }
    Node * node = inner->get_child(ind - 1);
    while (height > 1) {
        node = node->get_child(node->size());
        --height;
    }
    return dynamic_cast<Inner *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Inner * inner)
        -> Inner *
{
    if (!inner || !inner->m_parent) {
        return nullptr;
    }
    const Key & key = inner->get_key(0);
    Inner * parent = inner->m_parent;
    size_type height = 0;
    size_type ind = inner->size();
    while (parent && ind >= inner->size()) {
        inner = parent;
        ++height;
        parent = inner->m_parent;
        ind = inner->get_ind(key);
    }
    if (!parent && ind >= inner->size()) {
        return nullptr;
    }
    Node * node = inner->get_child(ind + 1);
    while (height > 1) {
        node = node->get_child(0);
        --height;
    }
    return dynamic_cast<Inner *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const_iterator begin, const_iterator end)
        -> iterator
{
    if (begin == end) {
        return iterator(end.m_current, end.m_pos);
    }
    std::vector<Key> m_keys;
    for (auto i = begin; i != end; ++i) {
        m_keys.push_back(i->first);
    }
    for (const Key & key : m_keys) {
        erase(key);
    }
    return upper_bound(m_keys.back());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const_iterator it)
        -> iterator
{
    if (it == end()) {
        return end();
    }
    Key key = it->first;
    auto [leaf, ind] = search(key);
    bool res = erase(leaf, ind, key);
    if (res) {
        --m_tree_size;
    }
    return upper_bound(key);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const Key & key)
        -> size_type
{
    auto [leaf, ind] = search(key);
    bool res = erase(leaf, ind, key);
    if (res) {
        --m_tree_size;
        return 1;
    }
    return 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::erase(Leaf * node, size_type ind, const Key & key)
{
    if (ind >= node->size() || less(key, node->get_key(ind))) {
        return false;
    }
    Leaf * left = find_left(node);
    Leaf * right = find_right(node);

    Key minimum = key;
    node->erase(ind, ind + 1);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        return true;
    }

    if (left && left->size() > m_order / 2 && !less(minimum, node->m_parent->get_key(0))) {
        node->shift(0);
        node->m_data[0] = left->m_data[left->size() - 1];
        left->erase(left->size() - 1, left->size());
        update(node, minimum, node->get_key(0));
    }
    else if (right && right->size() > m_order / 2) {
        node->shift(node->size());
        std::swap(node->m_data[node->size() - 1], right->m_data[0]);
        right->erase(0, 1);
        if (node->size() == 1) {
            update(node, minimum, node->get_key(0));
        }
        update(right, node->get_key(node->size() - 1), right->get_key(0));
    }
    else if (left && !less(minimum, node->m_parent->get_key(0))) {
        merge(left, node, minimum);
    }
    else if (right) {
        merge(node, right, right->get_key(0));
    }

    if (m_root->size() == 0 && m_root->get_child(0)) {
        m_root = m_root->get_child(0);
        m_root->m_parent->m_children[0] = nullptr;
        clear(m_root->m_parent);
        m_root->m_parent = nullptr;
    }

    return true;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::merge(Leaf * left, Leaf * right, const Key & key)
{
    for (size_type i = 0; i < right->size(); ++i) {
        left->shift(left->size());
        left->m_data[left->size() - 1] = right->m_data[i];
    }
    left->m_right = right->m_right;
    Inner * parent = right->m_parent;
    if (parent) {
        size_type ind = parent->get_ind(key);
        erase(parent, ind, key);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::erase(Inner * node, size_type ind, const Key & key)
{
    Inner * left = find_left(node);
    Inner * right = find_right(node);

    Key minimum = key;
    clear(node->get_child(ind));
    node->erase(ind - 1, ind);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        return;
    }

    if (left && left->size() > m_order / 2 && !less(minimum, node->m_parent->get_key(0))) {
        Leaf * node_min = find_min_leaf(node);

        node->shift(0);
        std::swap(node->m_children[0], node->m_children[1]);
        node->m_keys[0] = node_min->get_key(0);

        node->m_children[0] = left->m_children[left->size()];
        node->m_children[0]->m_parent = node;

        left->erase(left->size() - 1, left->size());
        update(node, minimum, find_min_leaf(node)->get_key(0));
    }
    else if (right && right->size() > m_order / 2) {
        node->shift(node->size());
        node->m_keys[node->size() - 1] = find_min_leaf(right)->get_key(0);
        node->m_children[node->size()] = right->m_children[0];
        node->m_children[node->size()]->m_parent = node;

        for (size_type i = 1; i < right->size(); ++i) {
            std::swap(right->m_keys[i - 1], right->m_keys[i]);
            std::swap(right->m_children[i - 1], right->m_children[i]);
        }
        std::swap(right->m_children[right->size() - 1], right->m_children[right->size()]);
        right->m_children[right->size()] = nullptr;
        --right->m_size;

        if (node->size() == 1) {
            update(node, minimum, find_min_leaf(node)->get_key(0));
        }
        update(right, node->m_keys[node->size() - 1], find_min_leaf(right)->get_key(0));
    }
    else if (left && !less(minimum, node->m_parent->get_key(0))) {
        merge(left, node, minimum);
    }
    else if (right) {
        merge(node, right, right->get_key(0));
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::merge(Inner * left, Inner * right, const Key & key)
{
    left->shift(left->size());
    left->m_keys[left->size() - 1] = find_min_leaf(right)->get_key(0);
    left->m_children[left->size()] = right->m_children[0];
    left->m_children[left->size()]->m_parent = left;
    right->m_children[0] = nullptr;

    for (size_type i = 0; i < right->size(); ++i) {
        left->shift(left->size());
        left->m_keys[left->size() - 1] = right->m_keys[i];
        left->m_children[left->size()] = right->m_children[i + 1];
        left->m_children[left->size()]->m_parent = left;
        right->m_children[i + 1] = nullptr;
    }

    Inner * parent = right->m_parent;
    if (parent) {
        size_type ind = parent->get_ind(key);
        erase(parent, ind, key);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::update(Node * node, const Key & prev, const Key & key)
{
    if (!node) {
        return;
    }

    while (node->m_parent) {
        size_type ind = node->m_parent->get_ind(prev);
        if (ind > 0) {
            node->m_parent->m_keys[ind - 1] = key;
            return;
        }
        node = node->m_parent;
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::clear()
{
    clear(m_root);
    m_root = new Leaf;
    m_tree_size = 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::clear(Node * node)
{
    if (node->get_child(0)) {
        for (size_type ind = 0; ind <= node->size(); ++ind) {
            clear(node->get_child(ind));
        }
    }
    delete node;
    node = nullptr;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::~BPTree()
{
    clear(m_root);
}
