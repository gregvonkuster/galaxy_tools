<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Samples for ${month_label}&nbsp;${year_label}
        </h3>
        <h4 align="center">
            Click Day to see samples for that day
        </h4>
        <table align="center" width="60%" class="colored">
            %if len( samples ) == 0:
                <tr>
                    <td colspan="2">
                        There are no samples for
                        ${month_label}&nbsp;${year_label}
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td class="half_width">Day</td>
                    <td class="half_width">Samples</td>
                </tr>
                <% ctr = 0 %>
                %for sample in samples:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>
                            <a href="${h.url_for( controller='samples', action='specified_date', specified_date=sample[0] )}">
                                ${sample[3]},
                                &nbsp;${month_label}&nbsp;${sample[1]},
                                &nbsp;${year_label}
                            </a>
                        </td>
                        <td>${sample[2]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
