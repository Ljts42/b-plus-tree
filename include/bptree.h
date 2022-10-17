#pragma once

#include <algorithm>
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

    virtual size_type size() const = 0;
    virtual const Key & get_key(size_type) const = 0;
    virtual size_type get_ind(const Key &) const = 0;
    virtual Node * get_child(size_type) const = 0;

    virtual ~Node() = default;
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Leaf : Node
{
    Leaf * m_right = nullptr;
    std::vector<value_type> m_data;

    Leaf()
        : m_data(0)
    {
        m_data.reserve(m_order);
    }

    Leaf(const Leaf & node)
        : m_data(node.m_data)
    {
    }

    size_type size() const override
    {
        return m_data.size();
    }

    const Key & get_key(size_type ind) const override
    {
        if (ind < size()) {
            return m_data[ind].first;
        }
        throw std::out_of_range("");
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
        return std::lower_bound(m_data.begin(), m_data.end(), key, compare) - m_data.begin();
    }

    ~Leaf()
    {
        this->m_parent = nullptr;
        m_right = nullptr;
    }
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Inner : Node
{
    std::vector<key_type> m_keys;
    std::vector<Node *> m_children;

    Inner()
        : m_keys(0)
        , m_children(0)
    {
        m_keys.reserve(m_order);
        m_children.reserve(m_order + 1);
    }

    Inner(const Inner & node)
        : m_keys(node.m_keys)
        , m_children(0)
    {
    }

    size_type size() const override
    {
        return m_keys.size();
    }

    const Key & get_key(size_type ind) const override
    {
        if (ind < size()) {
            return m_keys[ind];
        }
        throw std::out_of_range("");
    }

    Node * get_child(size_type ind) const override
    {
        return (ind < m_children.size()) ? m_children[ind] : nullptr;
    }

    size_type get_ind(const Key & key) const override
    {
        return std::upper_bound(m_keys.begin(), m_keys.end(), key, less) - m_keys.begin();
    }

    ~Inner()
    {
        this->m_parent = nullptr;
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
            inner->m_children.push_back(child);
            inner->m_children.back()->m_parent = inner;
        }
        return inner;
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::redirect()
{
    if (m_root && m_root->get_child(0)) {
        Leaf * leaf = find_min_leaf(m_root);
        if (!leaf->m_parent) {
            return;
        }
        Inner * node = leaf->m_parent;
        for (size_type ind = 0; ind < node->size(); ++ind) {
            leaf = dynamic_cast<Leaf *>(node->get_child(ind));
            leaf->m_right = dynamic_cast<Leaf *>(node->get_child(ind + 1));
        }
        Inner * next = find_right(node);
        while (next) {
            leaf = dynamic_cast<Leaf *>(node->get_child(node->size()));
            leaf->m_right = dynamic_cast<Leaf *>(next->get_child(0));
            node = next;
            next = find_right(node);
            for (size_type ind = 0; ind < node->size(); ++ind) {
                leaf = dynamic_cast<Leaf *>(node->get_child(ind));
                leaf->m_right = dynamic_cast<Leaf *>(node->get_child(ind + 1));
            }
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
    return (lower != end()) && !less(key, lower->first);
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
    if (ind >= leaf->size() || less(key, leaf->m_data[ind].first)) {
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
    if (ind >= leaf->size() || less(key, leaf->m_data[ind].first)) {
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
    if (ind >= leaf->size() || less(key, leaf->m_data[ind].first)) {
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
    leaf->m_data.insert(leaf->m_data.begin() + ind, std::make_pair(std::forward<K>(key), std::forward<V>(value)));
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
    std::move(left->m_data.begin() + m_order / 2, left->m_data.end(), std::back_inserter(right->m_data));
    left->m_data.erase(left->m_data.begin() + m_order / 2, left->m_data.end());
    right->m_right = left->m_right;
    left->m_right = right;

    if (!left->m_parent) {
        m_root = new Inner;
        Inner * parent = dynamic_cast<Inner *>(m_root);
        parent->m_children.push_back(left);
        left->m_parent = parent;
        parent->m_keys.push_back(right->get_key(0));
        parent->m_children.push_back(right);
        right->m_parent = parent;
    }
    else {
        Inner * parent = left->m_parent;
        right->m_parent = parent;
        size_type ind = std::lower_bound(parent->m_keys.begin(), parent->m_keys.end(), right->get_key(0), less) - parent->m_keys.begin();
        parent->m_keys.insert(parent->m_keys.begin() + ind, right->get_key(0));
        parent->m_children.insert(parent->m_children.begin() + ind + 1, right);
        if (parent->size() == m_order) {
            split(parent);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::split(Inner * left)
{
    Inner * right = new Inner;
    std::move(left->m_keys.begin() + m_order / 2 + 1, left->m_keys.end(), std::back_inserter(right->m_keys));
    left->m_keys.erase(left->m_keys.begin() + m_order / 2, left->m_keys.end());
    std::move(left->m_children.begin() + m_order / 2 + 1, left->m_children.end(), std::back_inserter(right->m_children));
    left->m_children.erase(left->m_children.begin() + m_order / 2 + 1, left->m_children.end());
    for (auto & i : right->m_children) {
        i->m_parent = right;
    }

    if (!left->m_parent) {
        m_root = new Inner;
        Inner * parent = dynamic_cast<Inner *>(m_root);
        parent->m_children.push_back(left);
        left->m_parent = parent;
        Leaf * rmin = find_min_leaf(right);
        parent->m_keys.push_back(rmin->get_key(0));
        parent->m_children.push_back(right);
        right->m_parent = parent;
    }
    else {
        Inner * parent = left->m_parent;
        right->m_parent = parent;
        Leaf * rmin = find_min_leaf(right);
        size_type ind = std::lower_bound(parent->m_keys.begin(), parent->m_keys.end(), rmin->get_key(0), less) - parent->m_keys.begin();
        parent->m_keys.insert(parent->m_keys.begin() + ind, rmin->get_key(0));
        parent->m_children.insert(parent->m_children.begin() + ind + 1, right);
        if (parent->size() == m_order) {
            split(parent);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Leaf * leaf)
        -> Leaf *
{
    if (!leaf || leaf->size() == 0 || !leaf->m_parent) {
        return nullptr;
    }
    const Key & key = leaf->get_key(0);
    Inner * parent = leaf->m_parent;
    Inner * node = parent;
    size_type height = 0;
    size_type ind = 0;
    while (parent && ind == 0) {
        node = parent;
        ++height;
        parent = node->m_parent;
        ind = node->get_ind(key);
    }
    if (!parent && ind == 0) {
        return nullptr;
    }
    if (height < 2) {
        return dynamic_cast<Leaf *>(node->get_child(ind - 1));
    }
    node = dynamic_cast<Inner *>(node->get_child(ind - 1));
    for (; height > 2; --height) {
        node = dynamic_cast<Inner *>(node->get_child(node->size()));
    }
    return dynamic_cast<Leaf *>(node->get_child(node->size()));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Inner * node)
        -> Inner *
{
    if (!node || node->size() == 0 || !node->m_parent) {
        return nullptr;
    }
    const Key & key = node->get_key(0);
    Inner * parent = node->m_parent;
    size_type height = 0;
    size_type ind = 0;
    while (parent && ind == 0) {
        node = parent;
        ++height;
        parent = node->m_parent;
        ind = node->get_ind(key);
    }
    if (!parent && ind == 0) {
        return nullptr;
    }
    node = dynamic_cast<Inner *>(node->get_child(ind - 1));
    for (; height > 1; --height) {
        node = dynamic_cast<Inner *>(node->get_child(node->size()));
    }
    return node;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Leaf * node)
        -> Leaf *
{
    return node->m_right;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Inner * node)
        -> Inner *
{
    if (!node || node->size() == 0 || !node->m_parent) {
        return nullptr;
    }
    const Key & key = node->get_key(0);
    Inner * parent = node->m_parent;
    size_type height = 0;
    size_type ind = node->size();
    while (parent && ind >= node->size()) {
        node = parent;
        ++height;
        parent = node->m_parent;
        ind = node->get_ind(key);
    }
    if (!parent && ind >= node->size()) {
        return nullptr;
    }
    node = dynamic_cast<Inner *>(node->get_child(ind + 1));
    for (; height > 1; --height) {
        node = dynamic_cast<Inner *>(node->get_child(0));
    }
    return node;
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
    node->m_data.erase(node->m_data.begin() + ind);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        if (ind == 0) {
            update(node, key, node->get_key(0));
        }
        return true;
    }

    if (left && left->size() > m_order / 2 && !less(minimum, node->m_parent->get_key(0))) {
        node->m_data.insert(node->m_data.begin(), left->m_data.back());
        left->m_data.pop_back();
        update(node, minimum, node->get_key(0));
    }
    else if (right && right->size() > m_order / 2) {
        node->m_data.insert(node->m_data.end(), right->m_data.front());
        right->m_data.erase(right->m_data.begin());
        if (node->size() == 1) {
            update(node, minimum, node->get_key(0));
        }
        update(right, node->m_data.back().first, right->get_key(0));
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
    std::move(right->m_data.begin(), right->m_data.end(), std::back_inserter(left->m_data));
    left->m_right = right->m_right;
    right->m_right = nullptr;
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
    node->m_keys.erase(node->m_keys.begin() + ind - 1);
    clear(node->get_child(ind));
    node->m_children.erase(node->m_children.begin() + ind);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        if (ind == 0) {
            update(node, key, find_min_leaf(node)->get_key(0));
        }
        return;
    }

    if (left && left->size() > m_order / 2 && !less(minimum, node->m_parent->get_key(0))) {
        node->m_keys.insert(node->m_keys.begin(), find_min_leaf(node)->get_key(0));
        left->m_keys.pop_back();
        node->m_children.insert(node->m_children.begin(), left->m_children.back());
        left->m_children.back()->m_parent = node;
        left->m_children.pop_back();
        update(node, minimum, find_min_leaf(node)->get_key(0));
    }
    else if (right && right->size() > m_order / 2) {
        node->m_keys.insert(node->m_keys.end(), find_min_leaf(right)->get_key(0));
        right->m_keys.erase(right->m_keys.begin());
        node->m_children.insert(node->m_children.end(), right->m_children.front());
        right->m_children.front()->m_parent = node;
        right->m_children.erase(right->m_children.begin());
        if (node->size() == 1) {
            update(node, minimum, find_min_leaf(node)->get_key(0));
        }
        update(right, node->m_keys.back(), find_min_leaf(right)->get_key(0));
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
    left->m_keys.push_back(find_min_leaf(right)->get_key(0));
    std::move(right->m_keys.begin(), right->m_keys.end(), std::back_inserter(left->m_keys));
    for (size_type i = 0; i <= right->size(); ++i) {
        left->m_children.push_back(right->get_child(i));
        left->m_children.back()->m_parent = left;
        right->m_children[i] = nullptr;
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
    Inner * parent = node->m_parent;
    if (!parent) {
        return;
    }
    while (parent) {
        size_type ind = parent->get_ind(prev);
        if (ind > 0) {
            parent->m_keys[ind - 1] = key;
            return;
        }
        node = parent;
        parent = node->m_parent;
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
