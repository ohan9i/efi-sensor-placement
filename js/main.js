import { initGame } from './game.js?v=1794145072';

function show(view) {
  const game = view !== 'learn';
  document.getElementById('gameView').hidden = !game;
  document.body.classList.toggle('dark', game);
  document.getElementById('top').hidden = game;
  document.querySelectorAll('.tabs button').forEach((b) => b.setAttribute('aria-selected', b.dataset.view === (game ? 'game' : 'learn')));
  window.scrollTo(0, 0);
  if (!game) window.dispatchEvent(new Event('learnshown'));
  else window.dispatchEvent(new Event('resize'));
}
window.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('.tabs button').forEach((b) => b.addEventListener('click', () => { history.replaceState(null, '', b.dataset.view === 'learn' ? '#learn' : '#'); show(b.dataset.view); }));
  initGame();
  show(location.hash === '#learn' ? 'learn' : 'game');
});
