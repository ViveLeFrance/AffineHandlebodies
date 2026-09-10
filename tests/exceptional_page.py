"""Exercise marker updates and stale replies in a JavaScript engine.

Run with python -m tests.exceptional_page. Requires the optional system
JavaScriptCore 4.1 introspection bindings; no browser or GUI is opened.
"""
import re
from pathlib import Path


def main():
    import gi
    gi.require_version('JavaScriptCore', '4.1')
    from gi.repository import JavaScriptCore as JSC
    ctx = JSC.Context.new()

    def run(code):
        result = ctx.evaluate(code, -1)
        error = ctx.get_exception()
        assert error is None, error.get_message()
        return result

    page = (Path(__file__).parents[1] / 'handlebodies/page.html').read_text()
    for script in re.findall(r'<script(?: [^>]*)?>(.*?)</script>', page, re.S):
        run('(function(){' + script + '\n})')
    run('''
      function check(value){if (!value) throw Error('assertion failed');}
      let body = {domain:'w-s**3-t*s', params:{}, pi:{expression:'t'}, rho:{expression:'w'}};
      function spec(){return body;}
      function $(id){return {addEventListener(){}};}
      const LAST = {pi:null}, EDITORS = {};
      const Plotly = {
        restyle(gd, update, indices){for (const i of indices)
          for (const k of Object.keys(update)) gd.data[i][k] = update[k][0];},
        relayout(gd, update){gd.relayout = update;}
      };
      const pending = [];
      function api(){return new Promise(resolve => pending.push(resolve));}
      function reply(result){pending.shift()({json:async () => result});}
      function cornerStyle(){return {};}
    ''')
    run(page[page.index('function equalScaleRanges('):page.index('function withOpeningView(')])
    run(page[page.index('const EXCEPTIONAL ='):page.index("$('fpi').addEventListener")])
    run(page[page.index('function editorTraces('):page.index('// The ends of a vanishing path')])
    run('''
      LAST.pi = {geometry: piGeometry(body), note: {}, gd: {isConnected:true,
        clientWidth:700, layout:{meta:{bounds:[-1,1,-1,1]}}, data:[exceptionalTrace([])]}};
      EDITORS.pi0 = {level:'pi', gd:{isConnected:true, layout:{},
        data:editorTraces([[0,0],[1,0]], [[1,0]], [0,0], 'pi')}};
      EDITORS.rho0 = {level:'rho', gd:{isConnected:true, layout:{},
        data:editorTraces([[0,0],[1,0]], [[1,0]], [0,0], 'rho')}};
      refreshExceptional();
      check(pending.length === 1);
      refreshExceptional();
      check(pending.length === 1);
      reply({status:'complete', values:[[3,4]]});
    ''')
    run('''
      check(EXCEPTIONAL.values.length === 1);
      check(LAST.pi.gd.data[0].text[0] === 'E0');
      check(EDITORS.pi0.gd.data[4].x[0] === 3);
      check(EDITORS.pi0.gd.data[0].text[0] === '0');
      check(EDITORS.rho0.gd.data.length === 4);
      check(LAST.pi.gd.relayout['xaxis.range'][1] >= 3);
      invalidateExceptional();
      check(EDITORS.pi0.gd.data[4].x.length === 0);
      refreshExceptional();
      invalidateExceptional();
      reply({status:'complete',values:[[99,99]]});
    ''')
    run('''
      check(EXCEPTIONAL.values.length === 0);
      refreshExceptional();
      reply({status:'timeout',values:[],message:'time limit'});
    ''')
    run('''
      check(LAST.pi.note.textContent.includes('not determined'));
      check(EXCEPTIONAL.values.length === 0);
      invalidateExceptional();
      body.pi.expression = 'other';
      refreshExceptional();
      check(pending.length === 0);
    ''')
    print('Page syntax, pi-only markers, labels, bounds, caching, stale replies and timeout display: ok')


if __name__ == '__main__':
    main()
