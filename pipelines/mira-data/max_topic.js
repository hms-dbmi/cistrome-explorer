a.forEach(d => { 
    let mi = 0;
    let mv = -1;
    for (let i = 0; i < c.length; i++) {
        if(mv < d[c[i]]) {
            mv = d[c[i]];
            mi = i;
        }
    }
    d['max_topic'] = c[mi]
});