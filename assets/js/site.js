// Topic filter on the publications page.
(function () {
  var bar = document.querySelector(".filter");
  if (!bar) return;
  bar.hidden = false;
  var buttons = bar.querySelectorAll("button");
  var items = document.querySelectorAll(".pubs .pub");
  var sections = document.querySelectorAll(".page > section");

  function apply(topic) {
    buttons.forEach(function (b) { b.setAttribute("aria-pressed", String(b.dataset.topic === topic)); });
    items.forEach(function (li) {
      var topics = (li.dataset.topics || "").split(" ");
      li.hidden = topic !== "all" && topics.indexOf(topic) === -1;
    });
    sections.forEach(function (s) {
      var visible = s.querySelectorAll(".pub:not([hidden])").length;
      var count = s.querySelector(".count");
      if (count) count.textContent = visible;
      s.hidden = visible === 0;
    });
    if (history.replaceState) history.replaceState(null, "", topic === "all" ? location.pathname : "#" + topic);
  }

  buttons.forEach(function (b) { b.addEventListener("click", function () { apply(b.dataset.topic); }); });
  var initial = location.hash.slice(1);
  if (initial && bar.querySelector('[data-topic="' + initial + '"]')) apply(initial);
})();
