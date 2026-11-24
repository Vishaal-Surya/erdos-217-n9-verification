import json, collections, sys
st = collections.Counter()
mism = []
sim = collections.defaultdict(set)
for l in open('verify/vb/ext_full.jsonl'):
    d = json.loads(l)
    st['factors'] += 1
    st['status:' + d['status']] += 1
    if d.get('nreal_vb') != d.get('nreal_builder'):
        mism.append(('nreal', d['key'], d.get('nreal_vb'), d.get('nreal_builder')))
    if d.get('pd_undecided'):
        mism.append(('pd_undecided', d['key']))
    if d['status'] == 'REALIZED':
        st['pd_embeddings_vb'] += d['pd_vb']
        st['pd_embeddings_builder'] += d['pd_builder']
        if d['pd_vb'] != d['pd_builder']:
            mism.append(('pd', d['key'], d['pd_vb'], d['pd_builder']))
        if not d.get('core_exact_ok'):
            mism.append(('core', d['key']))
        st['m=%d' % d['m']] += 1
    elif d.get('pd_vb'):
        mism.append(('nonrealized_with_pd_root', d['key'], d['status'], d['pd_vb']))
    for f in d.get('found') or []:
        st['found_numeric_sizes_ok'] += 1
        if f.get('exact') == 'REBUILD_FAILED':
            st['rebuild_failed'] += 1
            continue
        if f['sizes'] == list(range(1, 9)):
            st['exact_sizes_ok'] += 1
            if not f['collinear']:
                st['exact_sizes_ok_no_collinear'] += 1
                sig = tuple(round(x / max(f['dvals']), 4) for x in f['dvals'])
                sim[sig].add(d['key'])
            if f['crescent']:
                st['CRESCENT'] += 1
print(dict(st))
print('mismatches', len(mism), mism[:20])
print('no-collinear class-size-ok configs up to similarity:', len(sim))
for k, v in sim.items():
    print(k, sorted(v))
